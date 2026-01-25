"""
Jupiter magnetic field models

This module provides access to various Jupiter internal magnetic field models
based on Juno mission data and earlier observations.

# Available Models
- **JRM09**: Juno Reference Model through Perijove 9 (degree 10)
  - Reference: Connerney et al. (2018), Geophysical Research Letters
  - DOI: 10.1002/2018GL077312
  - Based on Juno's first 9 orbits

- **JRM33**: Juno Reference Model through Perijove 33 (degree up to 30)
  - Reference: Connerney et al. (2022), Journal of Geophysical Research: Planets
  - DOI: 10.1029/2021JE007055
  - Based on Juno's first 33 orbits
  - Well-determined through degree 13, useful information through degree 18

# Jupiter Physical Parameters
- Mean radius: 71,492 km (1 RJ)
- Reference frame: System III (1965) coordinates
- Dipole moment: ~4.17 Gauss RJ³
"""

# Jupiter physical constants
const JUPITER_RADIUS = 71492.0  # km (1 RJ)
const JUPITER_EPOCH = 2020.0    # Reference epoch (decimal year)

"""
    jupiter_model_info(model_name::String)

Get detailed information about a Jupiter magnetic field model.

# Arguments
- `model_name::String`: Model identifier

# Returns
- `Dict{String,Any}`: Model metadata including description, reference, DOI, degree, etc.

# Example
```julia
info = jupiter_model_info("jrm09")
println(info["description"])
println("Reference: ", info["reference"])
```
"""
function jupiter_model_info(model)
    model_name = lowercase(model)

    if model_name == "jrm09"
        return Dict{String, Any}(
            "name" => "JRM09",
            "description" => "Juno Reference Model through Perijove 9",
            "degree" => 10,
            "order" => 10,
            "reference" => "Connerney et al. (2018), Geophysical Research Letters",
            "doi" => "10.1002/2018GL077312",
            "year" => 2018,
            "mission_phase" => "First 9 Juno orbits",
            "notes" => "Well-validated standard model for Jupiter's internal field"
        )
    elseif model_name == "jrm33"
        return Dict{String, Any}(
            "name" => "JRM33",
            "description" => "Juno Reference Model through Perijove 33",
            "degree" => 30,
            "order" => 30,
            "recommended_degree" => 13,
            "reference" => "Connerney et al. (2022), Journal of Geophysical Research: Planets",
            "doi" => "10.1029/2021JE007055",
            "year" => 2022,
            "mission_phase" => "First 33 Juno orbits",
            "notes" => "Latest model. Coefficients well-determined through degree 13, " *
                "useful information through degree 18. Higher degrees less certain."
        )
    else
        error("Unknown model: $model_name")
    end
end
