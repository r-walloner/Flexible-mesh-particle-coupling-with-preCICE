#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/particle_handler.h>

#include <deal.II/particles/generators.h>

#include <deal.II/particles/data_out.h>

#include <precice/precice.hpp>

#include <cmath>
#include <iostream>
#include <vector>
#include <map>

namespace Step68
{
  using namespace dealii;

  enum class ParticipantRole
  {
    Fluid,
    Particle
  };

  ParticipantRole get_participant_role(const std::string &role)
  {
    if (role == "Fluid")
      return ParticipantRole::Fluid;
    else if (role == "Particle")
      return ParticipantRole::Particle;
    else
      AssertThrow(false, ExcMessage("Unknown participant role"));
  }

  class ParticleTrackingParameters : public ParameterAcceptor
  {
  public:
    ParticleTrackingParameters();

    std::string output_directory = "../solution/";

    unsigned int velocity_degree = 1;
    double time_step = 0.002;
    double final_time = 4.0;
    unsigned int output_interval = 10;
    unsigned int repartition_interval = 5;

    unsigned int fluid_refinement = 4;
    unsigned int particle_refinement = 4;
    unsigned int particle_insertion_refinement = 3;
  };

  ParticleTrackingParameters::ParticleTrackingParameters()
      : ParameterAcceptor("Particle Tracking Problem/")
  {
    add_parameter(
        "Velocity degree", velocity_degree, "", prm, Patterns::Integer(1));

    add_parameter("Output interval",
                  output_interval,
                  "Iteration interval between which output results are written",
                  prm,
                  Patterns::Integer(1));

    add_parameter("Repartition interval",
                  repartition_interval,
                  "Iteration interval at which the mesh is load balanced",
                  prm,
                  Patterns::Integer(1));

    add_parameter("Output directory", output_directory);

    add_parameter("Time step", time_step, "", prm, Patterns::Double());

    add_parameter("Final time",
                  final_time,
                  "End time of the simulation",
                  prm,
                  Patterns::Double());

    add_parameter("Fluid refinement",
                  fluid_refinement,
                  "Refinement level of the fluid domain",
                  prm,
                  Patterns::Integer(0));

    add_parameter("Particle grid refinement",
                  particle_refinement,
                  "Refinement level of the particle domain",
                  prm,
                  Patterns::Integer(0));

    add_parameter(
        "Particle insertion refinement",
        particle_insertion_refinement,
        "Refinement of the volumetric mesh used to insert the particles",
        prm,
        Patterns::Integer(0));
  }

  /* ------------------------------------------------------------------------
   * Fluid solver
   *
   *
   *
   *
   * ------------------------------------------------------------------------
   */

  template <int dim>
  class FluidSolver
  {
  public:
    FluidSolver(const ParticleTrackingParameters &par);
    void run();

  private:
    void setup_dofs();
    void setup_coupling();
    void interpolate_function_to_field();
    void output_field(const unsigned int it);
    void solve(const double dt);

    const ParticleTrackingParameters &par;

    MPI_Comm mpi_communicator;
    precice::Participant participant;
    ConditionalOStream pcout;

    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim> dh;
    const FESystem<dim> fe;
    MappingQ1<dim> mapping;

    LinearAlgebra::distributed::Vector<double> velocity_field;
    Functions::RayleighKotheVortex<dim> velocity;
  };

  template <int dim>
  FluidSolver<dim>::FluidSolver(const ParticleTrackingParameters &par)
      : par(par),
        mpi_communicator(MPI_COMM_WORLD),
        participant("Fluid",
                    "../precice-config.xml",
                    Utilities::MPI::this_mpi_process(mpi_communicator),
                    Utilities::MPI::n_mpi_processes(mpi_communicator)),
        pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
        triangulation(mpi_communicator),
        dh(triangulation),
        fe(FE_Q<dim>(par.velocity_degree) ^ dim),
        mapping(),
        velocity(4.0)
  {
  }

  template <int dim>
  void FluidSolver<dim>::setup_dofs()
  {
    GridGenerator::hyper_cube(triangulation, 0, 1);
    triangulation.refine_global(par.fluid_refinement);

    dh.distribute_dofs(fe);
    const IndexSet locally_owned_dofs = dh.locally_owned_dofs();
    const IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dh);

    velocity_field.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);
  }

  template <int dim>
  void FluidSolver<dim>::setup_coupling()
  {
    // Get the coordinates of the locally relevant support points.
    std::map<types::global_dof_index, Point<dim>> support_point_map{
        DoFTools::map_dofs_to_support_points(mapping, dh)};

    // Transform deal.II Points into format expected by preCICE
    std::vector<double> vertex_coordinates(support_point_map.size() * dim);
    unsigned int index = 0;
    for (const auto &point_pair : support_point_map)
    {
      const Point<dim> &point = point_pair.second;
      for (unsigned int d = 0; d < dim; ++d)
        vertex_coordinates[index++] = point[d];
    }

    // Give preCICE the vertex coordinates and initialize the participant
    std::vector<int> vertex_ids(dh.n_locally_owned_dofs());
    participant.setMeshVertices("Fluid-Mesh", vertex_coordinates, vertex_ids);

    participant.initialize();
  }

  template <int dim>
  void FluidSolver<dim>::interpolate_function_to_field()
  {
    velocity_field.zero_out_ghost_values();
    VectorTools::interpolate(mapping, dh, velocity, velocity_field);
    velocity_field.update_ghost_values();
  }

  template <int dim>
  void FluidSolver<dim>::output_field(const unsigned int it)
  {
    std::vector<std::string> solution_names(dim, "velocity");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
            dim, DataComponentInterpretation::component_is_part_of_vector);

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dh);
    data_out.add_data_vector(velocity_field,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches(mapping);

    const std::string output_folder(par.output_directory);
    const std::string file_name("fluid");

    pcout << "Writing fluid field file: " << file_name << '-' << it
          << std::endl;

    data_out.write_vtu_with_pvtu_record(
        output_folder, file_name, it, mpi_communicator, 6);
  }

  template <int dim>
  void FluidSolver<dim>::run()
  {
    DiscreteTime time(0, par.final_time, par.time_step);

    setup_dofs();
    setup_coupling();

    while (!time.is_at_end())
    {
      velocity.set_time(time.get_current_time());
      interpolate_function_to_field();

      if ((time.get_step_number() % par.output_interval) == 0)
        output_field(time.get_step_number());

      time.advance_time();
    }
  }

  /* ------------------------------------------------------------------------
   * Particle solver
   *
   *
   *
   *
   * ------------------------------------------------------------------------
   */

  template <int dim>
  class ParticleSolver
  {
  public:
    ParticleSolver(const ParticleTrackingParameters &par);
    void run();

  private:
    void generate_particles();
    void setup_coupling();
    void repartition();
    void euler_step(const double dt);
    void output_particles(const unsigned int it);

    unsigned int cell_weight(
        const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
        const CellStatus status) const;

    const ParticleTrackingParameters &par;

    MPI_Comm mpi_communicator;
    precice::Participant participant;
    ConditionalOStream pcout;

    parallel::distributed::Triangulation<dim> background_triangulation;
    Particles::ParticleHandler<dim> ph;
    MappingQ1<dim> mapping;
  };

  template <int dim>
  ParticleSolver<dim>::ParticleSolver(const ParticleTrackingParameters &par)
      : par(par),
        mpi_communicator(MPI_COMM_WORLD),
        participant("Particle",
                    "../precice-config.xml",
                    Utilities::MPI::this_mpi_process(mpi_communicator),
                    Utilities::MPI::n_mpi_processes(mpi_communicator)),
        pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
        background_triangulation(mpi_communicator)
  {
    background_triangulation.signals.weight.connect(
        [&](const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
            const CellStatus status) -> unsigned int
        {
          return this->cell_weight(cell, status);
        });
  }

  template <int dim>
  unsigned int ParticleSolver<dim>::cell_weight(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
      const CellStatus status) const
  {
    const unsigned int base_weight = 1;
    const unsigned int particle_weight = 10;

    unsigned int n_particles_in_cell = 0;
    switch (status)
    {
    case CellStatus::cell_will_persist:
    case CellStatus::cell_will_be_refined:
      n_particles_in_cell = ph.n_particles_in_cell(cell);
      break;

    case CellStatus::cell_invalid:
      break;

    case CellStatus::children_will_be_coarsened:
      for (const auto &child : cell->child_iterators())
        n_particles_in_cell += ph.n_particles_in_cell(child);
      break;

    default:
      DEAL_II_ASSERT_UNREACHABLE();
    }

    return base_weight + particle_weight * n_particles_in_cell;
  }

  template <int dim>
  void ParticleSolver<dim>::generate_particles()
  {
    ph.initialize(background_triangulation, mapping, 1 + dim);

    GridGenerator::hyper_cube(background_triangulation, 0, 1);
    background_triangulation.refine_global(par.particle_refinement);

    Point<dim> center;
    center[0] = 0.5;
    center[1] = 0.75;
    if (dim == 3)
      center[2] = 0.5;

    const double outer_radius = 0.15;
    const double inner_radius = 0.01;

    // Insert Particles using auxiliary triangulation
    parallel::distributed::Triangulation<dim> particle_triangulation(
        MPI_COMM_WORLD);
    GridGenerator::hyper_shell(
        particle_triangulation, center, inner_radius, outer_radius, 6);
    particle_triangulation.refine_global(par.particle_insertion_refinement);

    const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
        background_triangulation, IteratorFilters::LocallyOwnedCell());
    const auto global_bounding_boxes =
        Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

    std::vector<std::vector<double>> properties(
        particle_triangulation.n_locally_owned_active_cells(),
        std::vector<double>(dim + 1, 0.));

    Particles::Generators::quadrature_points(particle_triangulation,
                                             QMidpoint<dim>(),
                                             global_bounding_boxes,
                                             ph,
                                             mapping,
                                             properties);

    pcout << "Number of particles inserted: "
          << ph.n_global_particles() << std::endl;
  }

  template <int dim>
  void ParticleSolver<dim>::setup_coupling()
  {
    // TODO: fine-tune and verify this bounding box computation
    const BoundingBox<dim> bounding_box{
        GridTools::compute_mesh_predicate_bounding_box(
            background_triangulation, IteratorFilters::LocallyOwnedCell(), 2, true, 1)
            .front()};

    std::vector<double> bounding_box_vertices(dim * 2);
    for (unsigned int d = 0; d < dim; ++d)
    {
      bounding_box_vertices[d] = bounding_box.lower_bound(d);
      bounding_box_vertices[dim + d] = bounding_box.upper_bound(d);
    }

    participant.setMeshAccessRegion("Fluid-Mesh", bounding_box_vertices);

    participant.initialize();
  }

  template <int dim>
  void ParticleSolver<dim>::repartition()
  {
    ph.prepare_for_coarsening_and_refinement();
    background_triangulation.repartition();
    ph.unpack_after_coarsening_and_refinement();
  }

  template <int dim>
  void ParticleSolver<dim>::euler_step(const double dt)
  {
    const unsigned int this_mpi_rank = Utilities::MPI::this_mpi_process(mpi_communicator);
    Vector<double> particle_velocity(dim);

    for (auto &particle : ph)
    {
      Point<dim> &particle_location = particle.get_location();
      // TODO
      // particle_velocity = get from preCICE...
      particle_velocity = 0.0;

      for (int d = 0; d < dim; ++d)
        particle_location[d] += particle_velocity[d] * dt;

      ArrayView<double> properties = particle.get_properties();
      for (int d = 0; d < dim; ++d)
        properties[d] = particle_velocity[d];
      properties[dim] = this_mpi_rank;
    }
  }

  template <int dim>
  void ParticleSolver<dim>::output_particles(const unsigned int it)
  {
    Particles::DataOut<dim, dim> particle_output;

    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.emplace_back("process_id");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
            dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);

    particle_output.build_patches(ph,
                                  solution_names,
                                  data_component_interpretation);
    const std::string output_folder(par.output_directory);
    const std::string file_name("particles");

    pcout << "Writing particle output file: " << file_name << '-' << it
          << std::endl;

    particle_output.write_vtu_with_pvtu_record(
        output_folder, file_name, it, mpi_communicator, 6);
  }

  template <int dim>
  void ParticleSolver<dim>::run()
  {
    DiscreteTime time(0, par.final_time, par.time_step);

    generate_particles();
    setup_coupling();

    while (!time.is_at_end())
    {
      if ((time.get_step_number() % par.repartition_interval) == 0)
        repartition();

      euler_step(time.get_current_time());
      ph.sort_particles_into_subdomains_and_cells();

      if ((time.get_step_number() % par.output_interval) == 0)
        output_particles(time.get_step_number());

      time.advance_time();
    }
  }

} // namespace Step68

/* ------------------------------------------------------------------------
 * main
 *
 *
 *
 *
 * ------------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  using namespace Step68;
  using namespace dealii;
  deallog.depth_console(1);

  try
  {
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    // Read parameter file
    std::string prm_file;
    if (argc > 2)
      prm_file = argv[2];
    else
      prm_file = "../parameters.prm";
    ParticleTrackingParameters par;
    ParameterAcceptor::initialize(prm_file);

    // Determine the participant role and run
    if (argc <= 1)
    {
      std::cerr << "No participant role provided. Exiting..." << std::endl;
      return 1;
    }
    switch (get_participant_role(argv[1]))
    {
    case ParticipantRole::Fluid:
    {
      FluidSolver<2> fluid_solver(par);
      fluid_solver.run();
      break;
    }
    case ParticipantRole::Particle:
    {
      ParticleSolver<2> particle_solver(par);
      particle_solver.run();
      break;
    }
    }
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
