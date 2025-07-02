def generate(p):
    if p["solver"] == "AndersonJacksonFoam":
        return None
    else:
        return f"""
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvOptions;
}}

options
{{
    momentumSource
    {{
        type            vectorSemiImplicitSource;
        active          yes;
        
        vectorSemiImplicitSourceCoeffs
        {{
            selectionMode   all;
            volumeMode      specific;
            sources
            {{
                U
                {{
                    explicit
                    {{
                        type        exprField;
                        expression  "ExplicitMomentum{" / vol()" if p["write_mapping"] != "coarse-graining" else ""}";
                    }}

                    implicit
                    {{
                        type        exprField;
                        expression  "-ImplicitMomentum{" / vol()" if p["write_mapping"] != "coarse-graining" else ""}";
                    }}
                }}
            }}
        }}
    }}
}}
"""