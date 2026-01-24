"""
Coefficient loading and manipulation for MagneticModels.jl

This module provides functions to load Gauss coefficients from data files
and create GaussCoefficients objects.
"""

"""
    load_coefficients(filepath; max_degree = nothing)::GaussCoefficients

Load Gauss coefficients from a `filepath`.

- `max_degree::Union{Int,Nothing}`: Optional maximum degree to load (truncate higher degrees)

# Example
```julia
coeffs = load_coefficients("data/jupiter/jrm09.dat")
coeffs = load_coefficients("data/jupiter/jrm33.dat", max_degree=10)  # Truncate to degree 10
```
"""
function load_coefficients(filepath; max_degree = nothing)
    isfile(filepath) || error("File not found: $filepath")

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

    # Apply truncation if requested
    if max_degree !== nothing
        filter!(p -> p[1][1] <= max_degree, coeffs_dict)
    end

    # Build matrices (rest stays the same)
    degree = maximum(nm[1] for nm in keys(coeffs_dict))
    order = maximum(nm[2] for nm in keys(coeffs_dict))

    g_matrix = zeros(Float64, degree + 1, order + 1)
    h_matrix = zeros(Float64, degree + 1, order + 1)

    for ((n, m), (g, h)) in coeffs_dict
        g_matrix[n + 1, m + 1] = g
        h_matrix[n + 1, m + 1] = h
    end

    return GaussCoefficients(g_matrix, h_matrix, degree, order)
end

"""
    magnetic_dipole_moment(coeffs::GaussCoefficients, radius::Float64)

Calculate the magnetic dipole moment from g_1^0 coefficient.

The dipole moment is:
```math
M = -4π R³ g_1^0 / μ_0
```

where R is the planetary radius and μ_0 = 4π × 10^(-7) H/m.

# Arguments
- `coeffs::GaussCoefficients`: Gauss coefficients
- `radius::Float64`: Planetary radius in km

# Returns
- `Float64`: Magnetic dipole moment in units of 10^27 Am²

# Example
```julia
moment = magnetic_dipole_moment(coeffs, 71492.0)  # Jupiter
```
"""
function magnetic_dipole_moment(coeffs::GaussCoefficients, radius::Float64)
    g10 = coeffs.g[2, 1]  # g_1^0 (degree 1, order 0)

    # Convert radius from km to m
    R_m = radius * 1000.0

    # μ_0 = 4π × 10^(-7) H/m
    # M = -4π R³ g_1^0 / μ_0 = -R³ g_1^0 × 10^7
    # g_1^0 is in nT = 10^(-9) T
    # Result in Am² = 10^(-9) × R³ × 10^7 = R³ × 10^(-2)
    # Scale to 10^27 Am²: divide by 10^27

    # M in Am²
    M = -R_m^3 * g10 * 1.0e-2

    # Convert to 10^27 Am²
    return M / 1.0e27
end

"""
    truncate_coefficients(coeffs::GaussCoefficients, max_degree::Int)

Truncate coefficients to a maximum degree.

# Arguments
- `coeffs::GaussCoefficients`: Original coefficients
- `max_degree::Int`: Maximum degree to retain

# Returns
- `GaussCoefficients`: Truncated coefficients
"""
function truncate_coefficients(coeffs::GaussCoefficients, max_degree::Int)
    max_degree >= 1 || error("max_degree must be at least 1")
    max_degree >= coeffs.degree && return coeffs  # No truncation needed

    new_order = min(coeffs.order, max_degree)

    # Create new matrices
    g_new = coeffs.g[1:(max_degree + 1), 1:(new_order + 1)]
    h_new = coeffs.h[1:(max_degree + 1), 1:(new_order + 1)]

    return GaussCoefficients(g_new, h_new, max_degree, new_order)
end
