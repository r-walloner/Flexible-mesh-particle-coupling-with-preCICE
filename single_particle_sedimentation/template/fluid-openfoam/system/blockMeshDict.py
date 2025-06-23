def generate(p):
    return f"""
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale   1.0;

vertices
(
    ({-p["particle_diameter"] * 25 / 2} 0.00 {-p["particle_diameter"] * 25 / 2})
    ({p["particle_diameter"] * 25 / 2} 0.00 {-p["particle_diameter"] * 25 / 2})
    ({p["particle_diameter"] * 25 / 2} {p["particle_diameter"] * 75} {-p["particle_diameter"] * 25 / 2})
    ({-p["particle_diameter"] * 25 / 2} {p["particle_diameter"] * 75} {-p["particle_diameter"] * 25 / 2})
    ({-p["particle_diameter"] * 25 / 2} 0.00 {p["particle_diameter"] * 25 / 2})
    ({p["particle_diameter"] * 25 / 2} 0.00 {p["particle_diameter"] * 25 / 2})
    ({p["particle_diameter"] * 25 / 2} {p["particle_diameter"] * 75} {p["particle_diameter"] * 25 / 2})
    ({-p["particle_diameter"] * 25 / 2} {p["particle_diameter"] * 75} {p["particle_diameter"] * 25 / 2})
);

blocks
(
    hex (0 1 2 3 4 5 6 7)
    ({p["fluid_cells"][0]} {p["fluid_cells"][1]} {p["fluid_cells"][2]})
    simpleGrading (1 1 1)
);

boundary
(
    bottom
    {{
        type patch;
        faces
        (
            (0 1 5 4)
        );
    }}
    top
    {{
        type patch;
        faces
        (
            (2 3 7 6)
        );
    }}
    sides
    {{
        type wall;
        faces
        (
            (0 4 7 3)
            (1 2 6 5)
        );
    }}
    frontAndBack
    {{
        type wall;
        faces
        (
            (4 5 6 7)
            (0 3 2 1)
        );
    }}
);
"""