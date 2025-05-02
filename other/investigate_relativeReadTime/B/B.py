import precice

participant = precice.Participant("B", "../precice-config.xml", 0, 1)

participant.set_mesh_access_region("A-Mesh", [-1, 1, -1, 1])
participant.initialize()

t = 0
dt = 1

while participant.is_coupling_ongoing():
    t = t + dt
    print("Timestep: ", t)

    relative_read_time = dt * 0
    print("Relative read time: ", relative_read_time)
    print("Absolute read time: ", t - dt + relative_read_time)

    coordinates = [(0.0, 0.0)]
    data = participant.map_and_read_data("A-Mesh", "A-Data", coordinates, 0)
    print("Reading data: ", data)

    participant.advance(dt)

participant.finalize()