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
