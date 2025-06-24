def generate(p):
    return f"""
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      preciceDict;
}}

preciceConfig "../precice-config.xml";

participant Fluid;

modules (FF);

FF
{{
  nameAlpha alphaSolid;
  nameDragForce particleForce;
}};

interfaces
{{
  Interface1
  {{
    mesh              Fluid-Mesh;
    patches           ();
    locations         volumeCenters;

    readData ({"Alpha DragForce" if p["solver"] == "AndersonJacksonFoam" else "ExplicitMomentum ImplicitMomentum"});

    writeData
    (
      Velocity
    );
  }};
}};
"""