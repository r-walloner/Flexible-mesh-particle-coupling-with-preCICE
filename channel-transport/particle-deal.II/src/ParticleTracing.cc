#include <deal.II/base/array_view.h>
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include <precice/precice.hpp>

#include <iostream>
#include <string>
#include <vector>

namespace ParticleTracing
{
  using namespace dealii;
  using namespace precice;

  class Parameters : public ParameterAcceptor
  {
  public:
    Parameters();

    std::string precice_config_file = "../precice_config.xml";
    std::string precice_participant_name = "Particle";
    std::string precice_mesh_name = "Fluid-Mesh";
    std::string precice_data_name = "Velocity";

    std::string output_directory = "../solution/";
    std::string output_basename = "";

    int particle_grid_refinement = 4;

    double particle_insertion_x = 0.0;
    double particle_insertion_y = 0.0;
    double particle_insertion_z = 0.0;
    double particle_insertion_radius = 0.1;
    int particle_insertion_count = 1000;

    double start_time = 0.0;
    double final_time = 1.0;
    double time_step = 0.005;
    int output_interval = 10;
    int repartition_interval = 10;
  };

  Parameters::Parameters()
      : ParameterAcceptor("ParticleTracing")
  {
    add_parameter("precice_config_file", precice_config_file);
    add_parameter("precice_participant_name", precice_participant_name);
    add_parameter("precice_mesh_name", precice_mesh_name);
    add_parameter("precice_data_name", precice_data_name);
    add_parameter("output_directory", output_directory);
    add_parameter("output_basename", output_basename);
    add_parameter("particle_grid_refinement", particle_grid_refinement);
    add_parameter("particle_insertion_x", particle_insertion_x);
    add_parameter("particle_insertion_y", particle_insertion_y);
    add_parameter("particle_insertion_z", particle_insertion_z);
    add_parameter("particle_insertion_radius", particle_insertion_radius);
    add_parameter("particle_insertion_count", particle_insertion_count);
    add_parameter("start_time", start_time);
    add_parameter("final_time", final_time);
    add_parameter("time_step", time_step);
    add_parameter("output_interval", output_interval);
    add_parameter("repartition_interval", repartition_interval);
  }

  template <int dim>
  class InsertionDensityFunction : public Function<dim>
  {
  public:
    InsertionDensityFunction(const Point<dim> center, const double radius);
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

  private:
    const Point<dim> center = Point<dim>();
    const double radius = 0.0;
  };

  template <int dim>
  InsertionDensityFunction<dim>::InsertionDensityFunction(
      const Point<dim> center, const double radius)
      : Function<dim>(/*n_components=*/1),
        center(center),
        radius(radius)
  {
    AssertThrow(radius > 0.0, ExcMessage("Radius must be positive."));
  }

  template <int dim>
  double InsertionDensityFunction<dim>::value(
      const Point<dim> &p,
      const unsigned int /*component*/) const
  {
    double value = 0.0;
    if (p.distance(center) < radius)
      value = 1.0;
    return value;
  }

  template <int dim>
  class ParticleTracing
  {
  public:
    ParticleTracing(const Parameters &parameters);
    void run();

  private:
    void initialize();
    void insert_particles();
    void repartition();
    void step(const double dt);
    void output_particles(const unsigned int time_step) const;

    unsigned int cell_weight(
        const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
        const CellStatus &status) const;

    const Parameters &parameters;

    MPI_Comm mpi_comm;
    Participant precice;

    ConditionalOStream pcout;

    parallel::distributed::Triangulation<dim> grid;
    Particles::ParticleHandler<dim> particle_handler;
    MappingQ1<dim> mapping;
  };

  template <int dim>
  ParticleTracing<dim>::ParticleTracing(
      const Parameters &parameters)

      : parameters(parameters),
        mpi_comm(MPI_COMM_WORLD),
        precice(
            parameters.precice_participant_name,
            parameters.precice_config_file,
            Utilities::MPI::this_mpi_process(mpi_comm),
            Utilities::MPI::n_mpi_processes(mpi_comm)),
        pcout(
            std::cout,
            Utilities::MPI::this_mpi_process(mpi_comm) == 0),
        grid(mpi_comm)
  {
    grid.signals.weight.connect(
        [&](const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
            const CellStatus status) -> unsigned int
        {
          return this->cell_weight(cell, status);
        });
  }

  /**
   * @brief Compute the computational weight of a cell.
   *
   * The weight is used for load balancing during repartitioning of the grid.
   *
   * @tparam dim
   * @param cell The cell for which the weight is computed.
   * @param status The status of the cell.
   * @return The weight of the cell.
   *
   */
  template <int dim>
  unsigned int ParticleTracing<dim>::cell_weight(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
      const CellStatus &status) const
  {
    const unsigned int base_weight = 1;
    const unsigned int particle_weight = 10;

    unsigned int n_particles_in_cell = 0;
    switch (status)
    {
    case CellStatus::cell_will_persist:
    case CellStatus::cell_will_be_refined:
      n_particles_in_cell = particle_handler.n_particles_in_cell(cell);
      break;

    case CellStatus::cell_invalid:
      break;

    case CellStatus::children_will_be_coarsened:
      for (const auto &child : cell->child_iterators())
        n_particles_in_cell += particle_handler.n_particles_in_cell(child);
      break;

    default:
      DEAL_II_ASSERT_UNREACHABLE();
    }

    return base_weight + particle_weight * n_particles_in_cell;
  }

  /**
   * @brief Set up the particle handler, grid, and preCICE.
   *
   * @tparam dim
   */
  template <int dim>
  void ParticleTracing<dim>::initialize()
  {
    particle_handler.initialize(grid, mapping, dim + 1);

    // TODO: Read the grid boundary from a file and refine it
    GridGenerator::hyper_cube(grid, 0, 6);
    grid.refine_global(parameters.particle_grid_refinement);

    // Set up communication through preCICE
    // TODO: fine-tune and verify this bounding box computation
    const BoundingBox<dim> bounding_box{
        GridTools::compute_mesh_predicate_bounding_box(
            grid,
            IteratorFilters::LocallyOwnedCell(), 2, true, 1)
            .front()};

    // Convert the bounding box to the format expected by preCICE
    std::vector<double> bounding_box_vertices(dim * 2);
    for (unsigned int d = 0; d < dim; ++d)
    {
      bounding_box_vertices[2 * d] = bounding_box.lower_bound(d);
      bounding_box_vertices[2 * d + 1] = bounding_box.upper_bound(d);
    }

    // Set the bounding box for the participant and initialize
    const std::string mesh_name = parameters.precice_mesh_name;
    precice.setMeshAccessRegion(mesh_name, bounding_box_vertices);
    precice.initialize();

    // Print some information about the received mesh
    const int r_dim = precice.getMeshDimensions(mesh_name);
    const int r_vertex_count = precice.getMeshVertexSize(mesh_name);
    pcout << "Received mesh " << mesh_name << " has dimension " << r_dim
          << " and " << r_vertex_count << " vertices." << std::endl;
  }

  /**
   * @brief Create particles in the domain according to a probability density
   * function.
   *
   * @tparam dim
   */
  template <int dim>
  void ParticleTracing<dim>::insert_particles()
  {
    Point<dim> center;
    center[0] = parameters.particle_insertion_x;
    center[1] = parameters.particle_insertion_y;
    if (dim == 3)
      center[2] = parameters.particle_insertion_z;

    const double radius = parameters.particle_insertion_radius;

    // insert particles in using probability density function
    Particles::Generators::probabilistic_locations(
        grid,
        InsertionDensityFunction<dim>(center, radius),
        /*random_cell_selection=*/false,
        parameters.particle_insertion_count,
        particle_handler,
        mapping,
        /*random_number_seed=*/161);

    // set properties of each particle
    for (auto &particle : particle_handler)
    {
      ArrayView<double> properties = particle.get_properties();
      for (unsigned int p = 0; p < particle_handler.n_properties_per_particle(); ++p)
        properties[p] = 0.0;
    }

    pcout << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
  }

  /**
   * @brief Repartition the grid according to the cell_weight function.
   *
   * @tparam dim
   */
  template <int dim>
  void ParticleTracing<dim>::repartition()
  {
    pcout << "Repartitioning grid..." << std::endl;
    particle_handler.prepare_for_coarsening_and_refinement();
    grid.repartition();
    particle_handler.unpack_after_coarsening_and_refinement();
  }

  /**
   * @brief Compute the new position of each particle based on the fluid velocity
   * read from preCICE.
   *
   * @tparam dim
   * @param dt The time step size.
   *
   */
  template <int dim>
  void ParticleTracing<dim>::step(const double dt)
  {
    std::vector<double> location_vec(dim);
    std::vector<double> velocity(dim);
    const unsigned int this_mpi_rank = Utilities::MPI::this_mpi_process(mpi_comm);

    for (auto &particle : particle_handler)
    {
      ArrayView<double> properties = particle.get_properties();

      // get location
      Point<dim> &location = particle.get_location();
      for (int d = 0; d < dim; ++d)
        location_vec[d] = location[d];

      // get fluid velocity from preCICE
      precice.mapAndReadData(
          parameters.precice_mesh_name,
          parameters.precice_data_name,
          location_vec,
          /*relative_read_time=*/0,
          velocity);

      // update particle location and properties
      for (int d = 0; d < dim; ++d)
      {
        location[d] += velocity[d] * dt;
        properties[d] = velocity[d];
      }
      properties[3 * dim] = this_mpi_rank;

      particle_handler.sort_particles_into_subdomains_and_cells();
    }
  }

  /**
   * @brief Write the particles to file.
   *
   * @tparam dim
   * @param time_step The current time step. Used to name the output file.
   */
  template <int dim>
  void ParticleTracing<dim>::output_particles(const unsigned int time_step) const
  {
    Particles::DataOut<dim, dim> particle_output;

    // tell DataOut how to interpret the particle properties
    std::vector<std::string>
        data_component_name(particle_handler.n_properties_per_particle());
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(particle_handler.n_properties_per_particle());

    for (unsigned int p = 0; p < particle_handler.n_properties_per_particle(); ++p)
    {
      if (p < dim)
      {
        data_component_name[p] = "velocity";
        data_component_interpretation[p] =
            DataComponentInterpretation::component_is_part_of_vector;
      }
      else if (p == 3 * dim)
      {
        data_component_name[p] = "process_id";
        data_component_interpretation[p] =
            DataComponentInterpretation::component_is_scalar;
      }
    }

    particle_output.build_patches(particle_handler,
                                  data_component_name,
                                  data_component_interpretation);
    const std::string output_folder(parameters.output_directory);
    const std::string file_name = parameters.output_basename + "particles";

    pcout << "Writing particle output file: " << file_name << '-' << time_step
          << std::endl;

    particle_output.write_vtu_with_pvtu_record(
        output_folder, file_name, time_step, mpi_comm, 6);
  }

  /**
   * @brief Run the particle tracing.
   *
   * @tparam dim
   */
  template <int dim>
  void ParticleTracing<dim>::run()
  {
    DiscreteTime time(parameters.start_time, parameters.final_time);

    initialize();
    insert_particles();
    repartition();
    output_particles(0);

    while (precice.isCouplingOngoing())
    {
      time.set_next_step_size(
          std::min(parameters.time_step, precice.getMaxTimeStepSize()));
      pcout << "stepping dt = " << time.get_next_step_size() << std::endl;

      step(time.get_next_step_size());
      time.advance_time();

      if ((time.get_step_number() % parameters.repartition_interval) == 0)
        repartition();

      if ((time.get_step_number() % parameters.output_interval) == 0)
        output_particles(time.get_step_number());

      precice.advance(time.get_previous_step_size());
    }

    precice.finalize();
  }

}

int main(int argc, char *argv[])
{
  try
  {
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    // Read parameter file
    std::string parameter_file;
    if (argc > 1)
      parameter_file = argv[2];
    else
      parameter_file = "../parameters.prm";
    ParticleTracing::Parameters parameters;
    ParticleTracing::ParameterAcceptor::initialize(parameter_file);

    // Start the particle tracing
    ParticleTracing::ParticleTracing<2> particle_tracing(parameters);
    particle_tracing.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
