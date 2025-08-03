def generate(p):
    return f"""
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{{
    bottom
    {{
        type            fixedValue;
        value           uniform (0 {p["fluid_background_velocity"]} 0);
    }}
    spout
    {{
        type            fixedValue;
        value           uniform (0 {p["fluid_spout_velocity"]} 0);
    }}
    top
    {{
        type            zeroGradient;
    }}
    "sides|front|back"
    {{
        type            noSlip;
    }}
}}
"""