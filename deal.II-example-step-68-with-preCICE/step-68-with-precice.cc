#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
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

#include <algorithm>
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
    std::string output_basename = "";

    unsigned int velocity_degree = 1;
    double time_step = 0.002;
    double final_time = 4.0;
    unsigned int output_interval = 10;
    unsigned int repartition_interval = 5;

    unsigned int fluid_refinement = 4;
    unsigned int particle_refinement = 4;
    unsigned int particle_insertion_refinement = 3;

    std::string method = "euler_explicit";
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

    add_parameter("Output basename", output_basename);

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

    add_parameter("Method", method);
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
    FluidSolver(
        const ParticleTrackingParameters &par,
        std::string precice_config_file);
    void run();

  private:
    void setup_dofs();
    void setup_coupling();
    void solve();
    void output_field(const unsigned int it);
    void solve(const double dt);

    const ParticleTrackingParameters &par;

    MPI_Comm mpi_communicator;
    precice::Participant precice;
    ConditionalOStream pcout;

    double t;
    double dt;
    int step_number;

    parallel::distributed::Triangulation<dim> triangulation;
    std::vector<precice::VertexID> precice_vertex_ids;
    DoFHandler<dim> dh;
    const FESystem<dim> fe;
    MappingQ1<dim> mapping;

    LinearAlgebra::distributed::Vector<double> velocity_field;
    Functions::RayleighKotheVortex<dim> velocity;
  };

  template <int dim>
  FluidSolver<dim>::FluidSolver(
      const ParticleTrackingParameters &par,
      std::string precice_config_file)
      : par(par),
        mpi_communicator(MPI_COMM_WORLD),
        precice("Fluid",
                precice_config_file,
                Utilities::MPI::this_mpi_process(mpi_communicator),
                Utilities::MPI::n_mpi_processes(mpi_communicator)),
        pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
        t(0),
        dt(par.time_step),
        step_number(0),
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
    // Get the coordinates of the locally owned support points.
    const IndexSet local_dofs{dh.locally_owned_dofs()};
    const std::map<types::global_dof_index, Point<dim>> support_point_map{
        DoFTools::map_dofs_to_support_points(mapping, dh)};

    // Assemble data for preCICE coupling mesh
    std::vector<double> vertex_coordinates(local_dofs.size() / 2 * dim);
    unsigned int index = 0;
    auto dof_iterator = local_dofs.begin();
    while (dof_iterator != local_dofs.end())
    {
      for (unsigned int d = 0; d < dim; ++d)
        vertex_coordinates[index++] = support_point_map.at(*dof_iterator)[d];

      // we increment by 2 since there are 2 dofs at every vertex. This is verry hacky, but it works 😅
      dof_iterator++;
      dof_iterator++;
    }

    // Give preCICE the vertex coordinates and initialize the participant
    precice_vertex_ids.resize(vertex_coordinates.size() / dim);
    precice.setMeshVertices("Fluid-Mesh", vertex_coordinates, precice_vertex_ids);
    precice.initialize();
  }

  template <int dim>
  void FluidSolver<dim>::solve()
  {
    // Interpolate the solution to the analytical velocity function
    velocity_field.zero_out_ghost_values();
    VectorTools::interpolate(mapping, dh, velocity, velocity_field);
    velocity_field.update_ghost_values();

    // Evaluate the solution at the support points and send that data to preCICE
    std::vector<double> velocity_values(dh.n_locally_owned_dofs() / 2 * dim);

    double fluid_time = t;
    /*
      write(get_previous): initial, v_0, v_1
        -> works with read(dt) and get_previous
      write(get_current):  initial, v_1, v_2
        -> read(dt) and get_current: works
        -> read(dt) and get_previous: small error
        -> read(0) and get_curent: large error
        -> read(0) and get_previous: large error
      */

    velocity.set_time(fluid_time);
    pcout << "setting velocity function time to " << fluid_time
          << std::endl;
    std::vector<double> time_values(dh.n_locally_owned_dofs() / 2);
    std::fill(time_values.begin(), time_values.end(), fluid_time);

    const IndexSet local_dofs{dh.locally_owned_dofs()};
    const std::map<types::global_dof_index, Point<dim>> support_point_map{
        DoFTools::map_dofs_to_support_points(mapping, dh)};

    Functions::ConstantFunction<dim> test_velocity(fluid_time, dim);

    Vector<double> local_velocity(dim);
    unsigned int index = 0;
    auto dof_iterator = local_dofs.begin();
    while (dof_iterator != local_dofs.end())
    {
      velocity.vector_value(support_point_map.at(*dof_iterator), local_velocity);
      // test_velocity.vector_value(support_point_map.at(*dof_iterator), local_velocity);
      for (unsigned int d = 0; d < dim; ++d)
        velocity_values[index++] = local_velocity[d];

      // we increment by 2 since there are 2 dofs at every vertex. This is verry hacky, but it works 😅
      dof_iterator++;
      dof_iterator++;
    }

    precice.writeData("Fluid-Mesh", "Velocity", precice_vertex_ids, velocity_values);
    precice.writeData("Fluid-Mesh", "Time", precice_vertex_ids, time_values);
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
    const std::string file_name = par.output_basename + "fluid";

    pcout << "Writing fluid field file: " << file_name << '-' << it
          << std::endl;

    data_out.write_vtu_with_pvtu_record(
        output_folder, file_name, it, mpi_communicator, 6);
  }

  template <int dim>
  void FluidSolver<dim>::run()
  {
    setup_dofs();
    setup_coupling();
    output_field(0);

    while (precice.isCouplingOngoing())
    {
      dt = std::min(par.time_step, precice.getMaxTimeStepSize());
      t += dt;
      ++step_number;
      pcout << "step number " << step_number
            << "\n\tstepping " << dt
            << "\n\ttime is now at " << t
            << "\n\tdesired timestep was " << par.time_step
            << "\n\tprecice max timestep was " << precice.getMaxTimeStepSize()
            << std::endl;

      solve();

      if ((step_number % par.output_interval) == 0)
        output_field(step_number);

      precice.advance(dt);
    }

    precice.finalize();
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
    ParticleSolver(
        const ParticleTrackingParameters &par,
        std::string precice_config_file);
    void run();

  private:
    void generate_particles();
    void setup_coupling();
    void repartition();
    void step(const double dt);
    void step_euler_explicit(const double dt);
    void step_euler_implicit(const double dt);
    void step_trapezoidal(const double dt);
    void compute_error();
    void output_particles(const unsigned int it);

    unsigned int cell_weight(
        const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
        const CellStatus status) const;

    const ParticleTrackingParameters &par;

    MPI_Comm mpi_communicator;
    precice::Participant precice;
    ConditionalOStream pcout;

    double t;
    double dt;
    int step_number;

    parallel::distributed::Triangulation<dim> background_triangulation;
    Particles::ParticleHandler<dim> ph;
    MappingQ1<dim> mapping;

    // we need the same function that generates the fluid velocities in the fluid
    // "solver" here in the particle solver, to compute the exact particle locations
    // we compare against to quantify the positional error.
    Functions::RayleighKotheVortex<dim> fluid_velocity;
  };

  template <int dim>
  ParticleSolver<dim>::ParticleSolver(
      const ParticleTrackingParameters &par,
      std::string precice_config_file)
      : par(par),
        mpi_communicator(MPI_COMM_WORLD),
        precice("Particle",
                precice_config_file,
                Utilities::MPI::this_mpi_process(mpi_communicator),
                Utilities::MPI::n_mpi_processes(mpi_communicator)),
        pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
        t(0),
        dt(par.time_step),
        step_number(0),
        background_triangulation(mpi_communicator),
        fluid_velocity(4.0)
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
    ph.initialize(background_triangulation, mapping, 5 * dim + 1);

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
        std::vector<double>(3 * dim + 1, 0.));
    // first dim properties: interpolated velocity from preCICE
    // second dim properties: exact analytical position
    // third dim properties: exact analytical velocity
    // last 1 property: mpi rank

    Particles::Generators::quadrature_points(particle_triangulation,
                                             QMidpoint<dim>(),
                                             global_bounding_boxes,
                                             ph,
                                             mapping,
                                             properties);

    // set the initial analytical position of each particle
    for (auto &particle : ph)
    {
      const Point<dim> &location = particle.get_location();
      ArrayView<double> properties = particle.get_properties();
      for (int d = 0; d < dim; ++d)
        properties[d] = location[d];
    }

    pcout << "Number of particles inserted: "
          << ph.n_global_particles() << std::endl;
  }

  template <int dim>
  void ParticleSolver<dim>::setup_coupling()
  {
    // // TODO: fine-tune and verify this bounding box computation
    // const BoundingBox<dim> bounding_box{
    //     GridTools::compute_mesh_predicate_bounding_box(
    //         background_triangulation, IteratorFilters::LocallyOwnedCell(), 2, true, 1)
    //         .front()};

    // // Convert the bounding box to the format expected by preCICE
    // std::vector<double> bounding_box_vertices(dim * 2);
    // for (unsigned int d = 0; d < dim; ++d)
    // {
    //   bounding_box_vertices[2 * d] = bounding_box.lower_bound(d);
    //   bounding_box_vertices[2 * d + 1] = bounding_box.upper_bound(d);
    // }

    // Use the whole domain as the bounding box for now
    std::vector<double> bounding_box_vertices{0.0, 1.0, 0.0, 1.0};

    // Set the bounding box for the participant and initialize
    precice.setMeshAccessRegion("Fluid-Mesh", bounding_box_vertices);
    precice.initialize();

    // Print some information about the received mesh
    const int r_dim = precice.getMeshDimensions("Fluid-Mesh");
    const int r_vertex_count = precice.getMeshVertexSize("Fluid-Mesh");
    pcout << "Received mesh has dimension " << r_dim
          << " and " << r_vertex_count << " vertices." << std::endl;
  }

  template <int dim>
  void ParticleSolver<dim>::repartition()
  {
    ph.prepare_for_coarsening_and_refinement();
    background_triangulation.repartition();
    ph.unpack_after_coarsening_and_refinement();
  }

  template <int dim>
  void ParticleSolver<dim>::step(const double dt)
  {
    if (par.method == "euler_explicit")
      step_euler_explicit(dt);
    else if (par.method == "euler_implicit")
      step_euler_implicit(dt);
    else if (par.method == "trapezoidal")
      step_trapezoidal(dt);
    else
      AssertThrow(false, ExcMessage("Unknown method"));

    compute_error();

    ph.sort_particles_into_subdomains_and_cells();
  }

  template <int dim>
  void ParticleSolver<dim>::step_euler_explicit(const double dt)
  {
    pcout << "step_euler_explicit" << std::endl;

    std::vector<double> location_vec(dim);
    std::vector<double> velocity(dim);
    Point<dim> analytical_location{};
    Vector<double> analytical_velocity(dim);

    double fluid_time = t - dt;
    double relative_read_time = 0 * dt;

    fluid_velocity.set_time(fluid_time);
    pcout << "setting velocity function time to " << fluid_time
          << std::endl;
    pcout << "realtive read time is " << relative_read_time
          << std::endl;
    pcout << "absolute read time is " << (t - dt + relative_read_time)
          << std::endl;

    std::vector<double> read_time(1);
    std::vector<double> time_coordinates = {0.0, 0.0};
    precice.mapAndReadData(
        "Fluid-Mesh", "Time", time_coordinates, relative_read_time, read_time);
    pcout << "time read from the fluid participant is " << read_time.at(0)
          << std::endl;

    Functions::ConstantFunction<dim> test_fluid_velocity(fluid_time, dim);

    for (auto &particle : ph)
    {
      ArrayView<double> properties = particle.get_properties();

      // get location and analytical location
      Point<dim> &location = particle.get_location();
      for (int d = 0; d < dim; ++d)
      {
        location_vec[d] = location[d];
        analytical_location[d] = properties[d];
      }

      // get fluid velocity from preCICE and analytically
      precice.mapAndReadData(
          "Fluid-Mesh", "Velocity", location_vec, relative_read_time, velocity);
      fluid_velocity.vector_value(location, analytical_velocity);
      // test_fluid_velocity.vector_value(location, analytical_velocity);

      // update particle location and analytical location
      for (int d = 0; d < dim; ++d)
      {
        location[d] += velocity[d] * dt;
        analytical_location[d] += analytical_velocity[d] * dt;
      }

      // Update particle properties
      for (int d = 0; d < dim; ++d)
      {
        properties[d] = analytical_location[d];
        properties[2 * dim + d] = velocity[d];
        properties[3 * dim + d] = analytical_velocity[d];
      }
    }
  }

  template <int dim>
  void ParticleSolver<dim>::step_euler_implicit(const double dt)
  {
    pcout << "step_euler_implicit" << std::endl;

    std::vector<double> location_vec(dim);
    std::vector<double> velocity(dim);
    Point<dim> analytical_location{};
    Vector<double> analytical_velocity(dim);

    for (auto &particle : ph)
    {
      ArrayView<double> properties = particle.get_properties();

      // get location and analytical location
      Point<dim> &location = particle.get_location();
      for (int d = 0; d < dim; ++d)
      {
        location_vec[d] = location[d];
        analytical_location[d] = properties[d];
      }

      // get fluid velocity from preCICE and analytically
      precice.mapAndReadData(
          "Fluid-Mesh", "Velocity", location_vec, dt, velocity);
      fluid_velocity.set_time(t + dt);
      fluid_velocity.vector_value(location, analytical_velocity);

      // update particle location and analytical location
      for (int d = 0; d < dim; ++d)
      {
        location[d] += velocity[d] * dt;
        analytical_location[d] += analytical_velocity[d] * dt;
      }

      // Update particle properties
      for (int d = 0; d < dim; ++d)
      {
        properties[d] = analytical_location[d];
        properties[2 * dim + d] = velocity[d];
        properties[3 * dim + d] = analytical_velocity[d];
      }
    }
  }

  template <int dim>
  void ParticleSolver<dim>::step_trapezoidal(const double dt)
  {
    pcout << "step_trapezoidal" << std::endl;

    // We need to get velocities for all particles (first at t_n and then at t_n+1),
    // to utilize preCICE's caching. (see: https://precice.org/couple-your-code-just-in-time-mapping.html)
    std::map<types::particle_index, std::vector<double>> velocities_t0;
    std::map<types::particle_index, std::vector<double>> velocities_t1;

    for (Particles::ParticleAccessor<dim> &particle : ph)
    {
      const types::particle_index particle_id = particle.get_id();

      // get location
      std::vector<double> location_vec(dim);
      Point<dim> &location = particle.get_location();
      for (int d = 0; d < dim; ++d)
        location_vec[d] = location[d];

      // get fluid velocity from preCICE
      std::vector<double> velocity_t0(dim, 0.);
      std::vector<double> velocity_t1(dim, 0.);
      precice.mapAndReadData(
          "Fluid-Mesh", "Velocity", location_vec, 0, velocity_t0);
      precice.mapAndReadData(
          "Fluid-Mesh", "Velocity", location_vec, dt, velocity_t1);

      // store the velocities
      velocities_t0.insert({particle_id, velocity_t0});
      velocities_t1.insert({particle_id, velocity_t1});
    }

    // Now we can proceed to update each particle as usual
    Point<dim> analytical_location{};
    Vector<double> analytical_velocity_t0(dim);
    Vector<double> analytical_velocity_t1(dim);

    for (auto &particle : ph)
    {
      ArrayView<double> properties = particle.get_properties();
      types::particle_index particle_id = particle.get_id();

      // get analytical location
      Point<dim> &location = particle.get_location();
      for (int d = 0; d < dim; ++d)
        analytical_location[d] = properties[d];

      // get fluid velocity analytically
      fluid_velocity.set_time(t);
      fluid_velocity.vector_value(location, analytical_velocity_t0);
      fluid_velocity.set_time(t + dt);
      fluid_velocity.vector_value(location, analytical_velocity_t1);

      // update particle location and analytical location
      for (int d = 0; d < dim; ++d)
      {
        location[d] += (dt / 2) * (velocities_t0.at(particle_id)[d] +
                                   velocities_t1.at(particle_id)[d]);

        analytical_location[d] += (dt / 2) * (analytical_velocity_t0[d] +
                                              analytical_velocity_t1[d]);
      }

      
      // Update particle properties
      for (int d = 0; d < dim; ++d)
      {
        properties[d] = analytical_location[d];
        properties[2 * dim + d] = velocities_t1.at(particle_id)[d];
        properties[3 * dim + d] = analytical_velocity_t1[d];
      }
    }
  }

  template<int dim>
  void ParticleSolver<dim>::compute_error()
  {
    const unsigned int this_mpi_rank = Utilities::MPI::this_mpi_process(mpi_communicator);

    for (Particles::ParticleAccessor<dim> particle : ph)
    {
      Point<dim> &location = particle.get_location();
      ArrayView<double> properties = particle.get_properties();

      for (int d = 0; d < dim; ++d)
      {
        // see output_particle() for how properties are interpreted
        //
        // location_error   = location    - analytical_location
        properties[dim + d] = location[d] - properties[d];
        // velocity_error       = velocity                - analytical_velocity
        properties[4 * dim + d] = properties[2 * dim + d] - properties[3 * dim + d];
      }
      
      // process_id       = this_mpi_rank
      properties[5 * dim] = this_mpi_rank;
    }


  }

  template <int dim>
  void ParticleSolver<dim>::output_particles(const unsigned int it)
  {
    Particles::DataOut<dim, dim> particle_output;

    // tell DataOut how to interpret the particle properties
    std::vector<std::string>
        data_component_name(ph.n_properties_per_particle());
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(ph.n_properties_per_particle());

    for (unsigned int p = 0; p < ph.n_properties_per_particle(); ++p)
    {
      if (p < dim)
      {
        data_component_name[p] = "analytical_location";
        data_component_interpretation[p] =
        DataComponentInterpretation::component_is_part_of_vector;
      }
      else if (p < 2 * dim)
      {
        data_component_name[p] = "location_error";
        data_component_interpretation[p] =
            DataComponentInterpretation::component_is_part_of_vector;
      }
      else if (p < 3 * dim)
      {
        data_component_name[p] = "velocity";
        data_component_interpretation[p] =
            DataComponentInterpretation::component_is_part_of_vector;
      }
      else if (p < 4 * dim)
      {
        data_component_name[p] = "analytical_velocity";
        data_component_interpretation[p] =
            DataComponentInterpretation::component_is_part_of_vector;
      }
      else if (p < 5 * dim)
      {
        data_component_name[p] = "velocity_error";
        data_component_interpretation[p] =
            DataComponentInterpretation::component_is_part_of_vector;
      }
      else if (p == 5 * dim)
      {
        data_component_name[p] = "process_id";
        data_component_interpretation[p] =
            DataComponentInterpretation::component_is_scalar;
      }
    }

    particle_output.build_patches(ph,
                                  data_component_name,
                                  data_component_interpretation);
    const std::string output_folder(par.output_directory);
    const std::string file_name = par.output_basename + "particles";

    pcout << "Writing particle output file: " << file_name << '-' << it
          << std::endl;

    particle_output.write_vtu_with_pvtu_record(
        output_folder, file_name, it, mpi_communicator, 6);
  }

  template <int dim>
  void ParticleSolver<dim>::run()
  {
    generate_particles();
    repartition();
    output_particles(0);
    setup_coupling();

    while (precice.isCouplingOngoing())
    {
      dt = std::min(par.time_step, precice.getMaxTimeStepSize());
      t += dt;
      ++step_number;
      pcout << "step number " << step_number
            << "\n\tstepping " << dt
            << "\n\ttime is now at " << t
            << "\n\tdesired timestep was " << par.time_step
            << "\n\tprecice max timestep was " << precice.getMaxTimeStepSize()
            << std::endl;

      step(dt);

      if ((step_number % par.repartition_interval) == 0)
        repartition();

      if ((step_number % par.output_interval) == 0)
        output_particles(step_number);

      precice.advance(dt);
    }

    precice.finalize();
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

    // Read preCICE configuration file
    std::string precice_config_file;
    if (argc > 3)
      precice_config_file = argv[3];
    else
      precice_config_file = "../precice-config.xml";

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
      FluidSolver<2> fluid_solver(par, precice_config_file);
      fluid_solver.run();
      break;
    }
    case ParticipantRole::Particle:
    {
      ParticleSolver<2> particle_solver(par, precice_config_file);
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
