def generate(p):
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<precice-configuration experimental="True">

    <profiling mode="fundamental" />

    <data:vector name="Velocity" />
    <data:vector name="ExplicitMomentum" />
    <data:scalar name="ImplicitMomentum" />
    <data:scalar name="Alpha" />
    <data:vector name="DragForce" />

    <mesh name="Fluid-Mesh" dimensions="3">
        <use-data name="Velocity" />
        <use-data name="ExplicitMomentum" />
        <use-data name="ImplicitMomentum" />
        <use-data name="Alpha" />
        <use-data name="DragForce" />
    </mesh>

    <participant name="Fluid">
        <provide-mesh name="Fluid-Mesh" />
        <write-data name="Velocity" mesh="Fluid-Mesh" />
        <read-data name="ExplicitMomentum" mesh="Fluid-Mesh" />
        <read-data name="ImplicitMomentum" mesh="Fluid-Mesh" />
        <read-data name="Alpha" mesh="Fluid-Mesh" />
        <read-data name="DragForce" mesh="Fluid-Mesh" />
    </participant>

    <participant name="Particle">
        <receive-mesh name="Fluid-Mesh" from="Fluid" api-access="true" />
        <read-data name="Velocity" mesh="Fluid-Mesh" />
        <write-data name="ExplicitMomentum" mesh="Fluid-Mesh" />
        <write-data name="ImplicitMomentum" mesh="Fluid-Mesh" />
        <write-data name="Alpha" mesh="Fluid-Mesh" />
        <write-data name="DragForce" mesh="Fluid-Mesh" />

        {"""
        <mapping:nearest-neighbor direction="read" from="Fluid-Mesh" constraint="consistent" />
""" if p["read_mapping"] == "nearest-neighbor" else f"""
        <mapping:rbf direction="read" from="Fluid-Mesh" constraint="consistent">
            <basis-function:compact-polynomial-c4 support-radius="{p["read_mapping_radius"]}" />
        </mapping:rbf>
""" if p["read_mapping"] == "rbf" else ""}

        {"""
        <mapping:nearest-neighbor direction="write" to="Fluid-Mesh" constraint="conservative" />
""" if p["write_mapping"] == "nearest-neighbor" else f"""
        <mapping:coarse-graining
            direction="write"
            to="Fluid-Mesh"
            constraint="conservative"
            grain-dim="3"
            function-radius="{p["write_mapping_radius"]}" />
""" if p["write_mapping"] == "coarse-graining" else ""}
    </participant>

    <m2n:sockets acceptor="Fluid" connector="Particle" exchange-directory=".." />

    <coupling-scheme:serial-explicit>
        <participants first="Fluid" second="Particle" />
        <time-window-size value="{p["fluid_dt"]}" />
        <max-time value="{p["end_time"]}" />
        <exchange data="Velocity" mesh="Fluid-Mesh" from="Fluid" to="Particle" />
        <exchange data="ExplicitMomentum" mesh="Fluid-Mesh" from="Particle" to="Fluid" />
        <exchange data="ImplicitMomentum" mesh="Fluid-Mesh" from="Particle" to="Fluid" />
        <exchange data="Alpha" mesh="Fluid-Mesh" from="Particle" to="Fluid" />
        <exchange data="DragForce" mesh="Fluid-Mesh" from="Particle" to="Fluid" />
    </coupling-scheme:serial-explicit>

</precice-configuration>
"""