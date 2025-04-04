import precice

participant = precice.Participant("A", "../precice-config.xml", 0, 1)

mesh = [[0, 0]]
vertex_ids = participant.set_mesh_vertices("A-Mesh", mesh)
participant.initialize()

t = 0
dt = 1

while participant.is_coupling_ongoing():
    t = t + dt
    print("Timestep: ", t)

    data = t
    print("Writing data: ", data)
    participant.write_data("A-Mesh", "A-Data", vertex_ids, [data])

    participant.advance(dt)

participant.finalize()