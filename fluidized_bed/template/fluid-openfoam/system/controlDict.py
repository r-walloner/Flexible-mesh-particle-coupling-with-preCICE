def generate(p):
    return f"""
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     {p["solver"]};

startFrom       startTime;

startTime       0.0;

stopAt          endTime;

endTime         {p["end_time"]};

deltaT          {p["fluid_dt"]};

writeControl    adjustableRunTime;

writeInterval   {p["output_interval"]};

purgeWrite      0;

writePrecision  3;

writeCompression on;

timeFormat      general;

timePrecision   6;

libs ("libpreciceAdapterFunctionObject.so");
functions
{{
    preCICE_Adapter
    {{
        type preciceAdapterFunctionObject;
        errors strict;
    }}
}}
"""