import precice

participant = precice.Participant("B", "../precice-config.xml", 0, 1)

participant.set_mesh_access_region("A-Mesh", [-1, 1, -1, 1])

participant.initialize()

t = 0

while participant.is_coupling_ongoing():
    print("Timestep: ", t)
    print("Receiving data")
    data = participant.read_data("A-Mesh", "A-Data", [0], 0)
    print("Received data: ", t)

    t = t + 1
    participant.advance(1)

participant.finalize()