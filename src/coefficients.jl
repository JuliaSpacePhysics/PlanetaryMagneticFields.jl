"""
Coefficient loading and manipulation for PlanetaryMagneticFields.jl

This module provides functions to load Gauss coefficients from data files
and create GaussCoefficients objects.
"""

# Cache for loaded coefficients: filepath -> GaussCoefficients
const COEFFICIENT_CACHE = Dict{String, GaussCoefficients{Matrix{Float64}}}()

"""
    load_coefficients(filepath; use_cache=false)::GaussCoefficients

Load Gauss coefficients from the given `filepath`.

Set `use_cache=true` to cache loaded coefficients and reuse them on subsequent calls.

# Example
```julia
load_coefficients("data/jupiter/jrm09.dat")
```
"""
function load_coefficients(filepath; use_cache = false)
    isfile(filepath) || error("File not found: $filepath")

    # Check cache first
    if use_cache && haskey(COEFFICIENT_CACHE, filepath)
        return COEFFICIENT_CACHE[filepath]
    end

    # Storage: (n, m) -> (g, h)
    coeffs_dict = Dict{Tuple{Int, Int}, Tuple{Float64, Float64}}()

    open(filepath, "r") do file
        for line in eachline(file)
            (isempty(line) || startswith(line, '#')) && continue
            parts = split(line)
            coeff_type = parts[1]  # 'g' or 'h'
            n = parse(Int, parts[2])
            m = parse(Int, parts[3])
            value = parse(Float64, parts[4])

            # Get or initialize (g, h) pair
            current = get(coeffs_dict, (n, m), (0.0, 0.0))

            coeffs_dict[(n, m)] = if coeff_type == "g"
                (value, current[2])
            elseif coeff_type == "h"
                (current[1], value)
            end
        end
    end

    isempty(coeffs_dict) && error("No coefficients loaded")

    # Build matrices (rest stays the same)
    degree = maximum(nm[1] for nm in keys(coeffs_dict))
    order = maximum(nm[2] for nm in keys(coeffs_dict))

    g_matrix = zeros(Float64, degree + 1, order + 1)
    h_matrix = zeros(Float64, degree + 1, order + 1)

    for ((n, m), (g, h)) in coeffs_dict
        g_matrix[n + 1, m + 1] = g
        h_matrix[n + 1, m + 1] = h
    end

    result = GaussCoefficients(g_matrix, h_matrix)
    use_cache && (COEFFICIENT_CACHE[filepath] = result)
    return result
end

# Load all IGRF epoch files from a directory
function load_igrf_epochs(epochs = DateTime.(1900:5:2025))
    # Find all IGRF epoch files
    data_dir = pkgdir(@__MODULE__, "data/coeffs/earth")
    coefficients = map(epochs) do epoch
        filepath = joinpath(data_dir, "igrf$(year(epoch)).dat")
        load_coefficients(filepath)
    end
    return TimeVaryingGaussCoefficients(epochs, coefficients)
end
