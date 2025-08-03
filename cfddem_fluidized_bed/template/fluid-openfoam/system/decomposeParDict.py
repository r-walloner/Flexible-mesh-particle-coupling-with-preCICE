def generate(p):
    return f"""
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains {p["fluid_subdomains"][0] * p["fluid_subdomains"][1] * p["fluid_subdomains"][2]};

method simple;

simpleCoeffs
{{
    n ({p["fluid_subdomains"][0]} {p["fluid_subdomains"][1]} {p["fluid_subdomains"][2]});
}}

// ************************************************************************* //
"""