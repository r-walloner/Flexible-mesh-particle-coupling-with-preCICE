import precice

participant = precice.Participant("A", "../precice-config.xml", 0, 1)

mesh = [[0, 0]]
vertex_ids = participant.set_mesh_vertices("A-Mesh", mesh)

participant.initialize()

t = 0

while participant.is_coupling_ongoing():
    print("Timestep: ", t)
    print("Generating data")
    data = t
    print("Generated data: ", data)

    participant.write_data("A-Mesh", "A-Data", vertex_ids, [data])

    t = t + 1
    participant.advance(1)

participant.finalize()