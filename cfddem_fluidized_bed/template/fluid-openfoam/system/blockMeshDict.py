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
    (-0.075 0.00 -0.0075)
    (-0.005 0.00 -0.0075)
    (-0.005 0.75 -0.0075)
    (-0.075 0.75 -0.0075)
    (-0.075 0.00 +0.0075)
    (-0.005 0.00 +0.0075)
    (-0.005 0.75 +0.0075)
    (-0.075 0.75 +0.0075)

    (-0.005 0.00 -0.0075)
    (+0.005 0.00 -0.0075)
    (+0.005 0.75 -0.0075)
    (-0.005 0.75 -0.0075)
    (-0.005 0.00 +0.0075)
    (+0.005 0.00 +0.0075)
    (+0.005 0.75 +0.0075)
    (-0.005 0.75 +0.0075)

    (+0.005 0.00 -0.0075)
    (+0.075 0.00 -0.0075)
    (+0.075 0.75 -0.0075)
    (+0.005 0.75 -0.0075)
    (+0.005 0.00 +0.0075)
    (+0.075 0.00 +0.0075)
    (+0.075 0.75 +0.0075)
    (+0.005 0.75 +0.0075)
);

blocks
(
    hex (0 1 2 3 4 5 6 7)
    ({int(p["fluid_cells"][0] * (70 / 150))} {p["fluid_cells"][1]} {p["fluid_cells"][2]})
    simpleGrading (1 1 1)

    hex (8 9 10 11 12 13 14 15)
    ({int(p["fluid_cells"][0] * (10 / 150))} {p["fluid_cells"][1]} {p["fluid_cells"][2]})
    simpleGrading (1 1 1)

    hex (16 17 18 19 20 21 22 23)
    ({int(p["fluid_cells"][0] * (70 / 150))} {p["fluid_cells"][1]} {p["fluid_cells"][2]})
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
            (16 17 21 20)
        );
    }}
    spout
    {{
        type patch;
        faces
        (
            (8 9 13 12)
        );
    }}
    top
    {{
        type patch;
        faces
        (
            (2 3 7 6)
            (10 11 15 14)
            (18 19 23 22)
        );
    }}
    sides
    {{
        type wall;
        faces
        (
            (0 4 7 3)
            (17 18 22 21)
        );
    }}
    front
    {{
        type wall;
        faces
        (
            (4 5 6 7)
            (12 13 14 15)
            (20 21 22 23)
        );
    }}
    back
    {{
        type wall;
        faces
        (
            (0 3 2 1)
            (8 11 10 9)
            (16 19 18 17)
        );
    }}

    // patches to merge blocks on
    merge0l
    {{
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }}
    merge0r
    {{
        type patch;
        faces
        (
            (8 12 15 11)
        );
    }}
    merge1l
    {{
        type patch;
        faces
        (
            (9 10 14 13)
        );
    }}
    merge1r
    {{
        type patch;
        faces
        (
            (16 20 23 19)
        );
    }}
);

mergePatchPairs
(
    (merge0l merge0r)
    (merge1l merge1r)
);
"""