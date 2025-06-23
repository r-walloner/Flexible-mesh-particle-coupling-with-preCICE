def generate(p):
    return f"""
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}}

transportModel  Newtonian;

nu              {p["fluid_viscosity"] / p["fluid_density"]}; // Ns/m^3
rho             {p["fluid_density"]}; // kg/m^3
g               (0 -9.81 0);
"""
