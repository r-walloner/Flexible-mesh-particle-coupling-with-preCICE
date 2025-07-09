def generate(p):
    if p["read_mapping"] == "nearest-neighbor":
        read_mapping = """
        <mapping:nearest-neighbor direction="read" from="Fluid-Mesh" constraint="consistent" />
"""
    
    elif p["read_mapping"] == "rbf":
        read_mapping = f"""
        <mapping:rbf direction="read" from="Fluid-Mesh" constraint="consistent">
            <basis-function:compact-polynomial-c4 support-radius="{p["read_mapping_radius"]}" />
        </mapping:rbf>
"""
        
    if p["write_mapping"] == "nearest-neighbor":
        write_mapping = """
        <mapping:nearest-neighbor direction="write" to="Fluid-Mesh" constraint="conservative" />
"""

    elif p["write_mapping"] == "rbf":
        write_mapping = f"""
        <mapping:rbf direction="write" to="Fluid-Mesh" constraint="conservative">
            <basis-function:compact-polynomial-c6 support-radius="{p["write_mapping_radius"]}"/>
        </mapping:rbf>
"""
        
    elif p["write_mapping"] == "coarse-graining":
        write_mapping = f"""
        <mapping:coarse-graining
            direction="write"
            to="Fluid-Mesh"
            constraint="conservative"
            grain-dim="3"
            function-radius="{p["write_mapping_radius"]}" />
"""
        
    if p["coupling"] == "explicit":
        data = """
    <data:vector name="DragForce" />
"""

        mesh_use_data = """
        <use-data name="DragForce" />
"""

        fluid_read_data = """
        <read-data name="DragForce" mesh="Fluid-Mesh" />
"""

        particle_write_data = """
        <write-data name="DragForce" mesh="Fluid-Mesh" />
"""

        exchange_data = """
        <exchange data="DragForce" mesh="Fluid-Mesh" from="Particle" to="Fluid" />
"""

    

    elif p["coupling"] == "semi_implicit":
        data = """
    <data:scalar name="ImplicitMomentum" />
    <data:vector name="ExplicitMomentum" />
"""

        mesh_use_data = """
        <use-data name="ImplicitMomentum" />
        <use-data name="ExplicitMomentum" />
"""

        fluid_read_data = """
        <read-data name="ImplicitMomentum" mesh="Fluid-Mesh" />
        <read-data name="ExplicitMomentum" mesh="Fluid-Mesh" />
"""

        particle_write_data = """
        <write-data name="ImplicitMomentum" mesh="Fluid-Mesh" />
        <write-data name="ExplicitMomentum" mesh="Fluid-Mesh" />
"""

        exchange_data = """
        <exchange data="ImplicitMomentum" mesh="Fluid-Mesh" from="Particle" to="Fluid" />
        <exchange data="ExplicitMomentum" mesh="Fluid-Mesh" from="Particle" to="Fluid" />
"""

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<precice-configuration experimental="True">

    <log>
        <sink
            type="stream"
            output="stdout"
            filter="{"%Severity% >= debug" if p["precice_debug_log"] else "(%Severity% > debug) and not ((%Severity% = info) and (%Rank% != 0))"}" />
    </log>

    <profiling mode="{p["precice_profiling"]}" />

    <data:vector name="Velocity" />
    <data:scalar name="Alpha" />
    {data}

    <mesh name="Fluid-Mesh" dimensions="3">
        <use-data name="Velocity" />
        <use-data name="Alpha" />
        {mesh_use_data}
    </mesh>

    <participant name="Fluid">
        <provide-mesh name="Fluid-Mesh" />
        <write-data name="Velocity" mesh="Fluid-Mesh" />
        <read-data name="Alpha" mesh="Fluid-Mesh" />
        {fluid_read_data}
    </participant>

    <participant name="Particle">
        <receive-mesh name="Fluid-Mesh" from="Fluid" api-access="true" />
        <read-data name="Velocity" mesh="Fluid-Mesh" />
        <write-data name="Alpha" mesh="Fluid-Mesh" />
        {particle_write_data}
        {read_mapping}
        {write_mapping}
    </participant>

    <m2n:sockets acceptor="Fluid" connector="Particle" exchange-directory=".." />

    <coupling-scheme:serial-explicit>
        <participants first="Fluid" second="Particle" />
        <time-window-size value="{p["fluid_dt"]}" />
        <max-time value="{p["end_time"]}" />
        <exchange data="Velocity" mesh="Fluid-Mesh" from="Fluid" to="Particle" />
        <exchange data="Alpha" mesh="Fluid-Mesh" from="Particle" to="Fluid" />
        {exchange_data}
    </coupling-scheme:serial-explicit>

</precice-configuration>
"""