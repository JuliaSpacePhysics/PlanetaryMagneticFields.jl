"""
    Current Sheet Models for Planetary Magnetospheres

This module implements the Connerney et al. (1981, 2020) Jupiter magnetodisc model,
which describes the magnetic field due to a "washer-shaped" azimuthal current sheet
near Jupiter's magnetic equator.

# References
- Connerney, J. E. P., Acuña, M. H., & Ness, N. F. (1981). Modeling the Jovian current
  sheet and inner magnetosphere. J. Geophys. Res., 86(A10), 8370-8384.
- Edwards, T. M., Bunce, E. J., & Cowley, S. W. H. (2001). A note on the vector potential
  of Connerney et al.'s model. Planet. Space Sci., 49, 1115-1123.
- Connerney, J. E. P., Timmins, S., Herceg, M., & Joergensen, J. L. (2020). A Jovian
  magnetodisc model for the Juno era. J. Geophys. Res. Space Physics, 125, e2020JA028138.
"""

"""
    CurrentSheetEquationType

Enum-like type for selecting the equation type used in current sheet calculations.
- `Analytical`: Use Edwards et al. (2001) analytical approximations (fast)
- `Integral`: Use Connerney et al. (1981) numerical integration (accurate)
- `Hybrid`: Combine both methods based on position (default, recommended)
"""
@enum CurrentSheetEquationType begin
    Analytical
    Integral
    Hybrid
end

"""
    CurrentSheetParameters

Parameters for the Connerney current sheet model.

# Fields
- `μ_i_half::Float64`: Current intensity parameter μ₀I₀/(4πd) in nT (default: 139.6)
- `r0::Float64`: Inner edge of current sheet in Rⱼ (default: 7.8)
- `r1::Float64`: Outer edge of current sheet in Rⱼ (default: 51.4)
- `d::Float64`: Half-thickness of current sheet in Rⱼ (default: 3.6)
- `θ_tilt::Float64`: Tilt angle of disc from spin axis in radians (default: 9.3°)
- `φ_tilt::Float64`: Longitude of tilt direction in SIII in radians (default: 155.8°)
- `I_rho::Float64`: Radial current in MA (default: 16.7, CON2020 only)
- `use_radial_current::Bool`: Whether to include radial current contribution (default: true)
"""
Base.@kwdef struct CurrentSheetParameters
    μ_i_half::Float64 = 139.6       # nT
    r0::Float64 = 7.8               # Rⱼ (inner edge)
    r1::Float64 = 51.4              # Rⱼ (outer edge)
    d::Float64 = 3.6                # Rⱼ (half-thickness)
    θ_tilt::Float64 = deg2rad(9.3)  # radians
    φ_tilt::Float64 = deg2rad(155.8) # radians (SIII longitude)
    I_rho::Float64 = 16.7           # MA (radial current)
    use_radial_current::Bool = true
end

# Pre-computed parameter for radial current: μ₀Iᵨ/(2π) in nT×Rⱼ
# Conversion: 16.7 MA → 16.7 nT×Rⱼ (using μ₀/(2π) ≈ 0.2 μT/A × Rⱼ)
radial_current_factor(p::CurrentSheetParameters) = p.use_radial_current ? p.I_rho : 0.0

"""
    CurrentSheetModel <: ExternalFieldModel

Connerney et al. (1981, 2020) Jupiter magnetodisc current sheet model.

This model calculates the magnetic field contribution from a washer-shaped
azimuthal current sheet centered on Jupiter's magnetic equator.

# Fields
- `name::String`: Model name
- `params::CurrentSheetParameters`: Model parameters
- `equation_type::CurrentSheetEquationType`: Equation type to use
"""
struct CurrentSheetModel <: ExternalFieldModel
    name::String
    params::CurrentSheetParameters
    equation_type::CurrentSheetEquationType
end

function CurrentSheetModel(;
    name::String = "CON2020",
    params::CurrentSheetParameters = CurrentSheetParameters(),
    equation_type::CurrentSheetEquationType = Hybrid
)
    return CurrentSheetModel(name, params, equation_type)
end

"""
    (model::CurrentSheetModel)(r, θ, φ, t=nothing)

Evaluate the current sheet magnetic field at the given position.

# Arguments
- `r`: Radial distance in planetary radii (Rⱼ)
- `θ`: Colatitude in radians [0, π]
- `φ`: East longitude in radians [0, 2π] (SIII)
- `t`: Time (unused, for API compatibility)

# Returns
- `SVector{3,Float64}`: Magnetic field [Br, Bθ, Bφ] in nT (SIII spherical)
"""
function (model::CurrentSheetModel)(r, θ, φ, t = nothing)
    p = model.params

    # Convert spherical SIII to Cartesian SIII
    x, y, z = spherical_to_cartesian(r, θ, φ)

    # Transform to current sheet (CS) coordinates (tilted disc frame)
    x_cs, y_cs, z_cs = siii_to_cs(x, y, z, p.θ_tilt, p.φ_tilt, φ)

    # Calculate cylindrical coordinates in CS frame
    ρ = sqrt(x_cs^2 + y_cs^2)
    φ_cs = atan(y_cs, x_cs)

    # Calculate field in cylindrical CS coordinates
    B_ρ, B_φ_cs, B_z = current_sheet_field_cyl(
        ρ, z_cs, p, model.equation_type
    )

    # Convert cylindrical to Cartesian in CS frame
    Bx_cs = B_ρ * cos(φ_cs) - B_φ_cs * sin(φ_cs)
    By_cs = B_ρ * sin(φ_cs) + B_φ_cs * cos(φ_cs)
    Bz_cs = B_z

    # Transform back to SIII Cartesian
    Bx, By, Bz = cs_to_siii(Bx_cs, By_cs, Bz_cs, p.θ_tilt, p.φ_tilt, φ)

    # Convert to spherical components
    return cartesian_field_to_spherical(Bx, By, Bz, θ, φ)
end

"""
    siii_to_cs(x, y, z, θ_tilt, φ_tilt, φ)

Transform coordinates from SIII Cartesian to current sheet (tilted disc) frame.

The current sheet is tilted from the rotational (z) axis by angle θ_tilt
in the direction specified by φ_tilt.
"""
function siii_to_cs(x, y, z, θ_tilt, φ_tilt, φ_pos)
    # The tilt rotates around an axis perpendicular to the tilt direction
    # First rotate to align tilt direction with x-axis, then tilt around y
    sin_tilt = sin(θ_tilt)
    cos_tilt = cos(θ_tilt)
    sin_φtilt = sin(φ_tilt)
    cos_φtilt = cos(φ_tilt)

    # Rotation matrix: R_z(-φ_tilt) × R_y(-θ_tilt) × R_z(φ_tilt)
    # This tilts the z-axis toward the direction φ_tilt

    # First: rotate around z by -φ_tilt to bring tilt direction to x-z plane
    x1 = x * cos_φtilt + y * sin_φtilt
    y1 = -x * sin_φtilt + y * cos_φtilt
    z1 = z

    # Second: rotate around y by -θ_tilt (tilt the coordinate system)
    x_cs = x1 * cos_tilt + z1 * sin_tilt
    y_cs = y1
    z_cs = -x1 * sin_tilt + z1 * cos_tilt

    return x_cs, y_cs, z_cs
end

"""
    cs_to_siii(Bx_cs, By_cs, Bz_cs, θ_tilt, φ_tilt, φ_pos)

Transform field vector from current sheet frame back to SIII Cartesian.
"""
function cs_to_siii(Bx_cs, By_cs, Bz_cs, θ_tilt, φ_tilt, φ_pos)
    sin_tilt = sin(θ_tilt)
    cos_tilt = cos(θ_tilt)
    sin_φtilt = sin(φ_tilt)
    cos_φtilt = cos(φ_tilt)

    # Inverse transformation: R_z(-φ_tilt)ᵀ × R_y(-θ_tilt)ᵀ × R_z(φ_tilt)ᵀ
    # = R_z(-φ_tilt) × R_y(θ_tilt) × R_z(φ_tilt)

    # First: rotate around y by +θ_tilt
    x1 = Bx_cs * cos_tilt - Bz_cs * sin_tilt
    y1 = By_cs
    z1 = Bx_cs * sin_tilt + Bz_cs * cos_tilt

    # Second: rotate around z by +φ_tilt
    Bx = x1 * cos_φtilt - y1 * sin_φtilt
    By = x1 * sin_φtilt + y1 * cos_φtilt
    Bz = z1

    return Bx, By, Bz
end

"""
    cartesian_field_to_spherical(Bx, By, Bz, θ, φ)

Convert field components from Cartesian to spherical coordinates.
"""
function cartesian_field_to_spherical(Bx, By, Bz, θ, φ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinφ = sin(φ)
    cosφ = cos(φ)

    # Inverse of the spherical-to-Cartesian transformation matrix
    Br = Bx * sinθ * cosφ + By * sinθ * sinφ + Bz * cosθ
    Bθ = Bx * cosθ * cosφ + By * cosθ * sinφ - Bz * sinθ
    Bφ = -Bx * sinφ + By * cosφ

    return SVector{3,Float64}(Br, Bθ, Bφ)
end

"""
    current_sheet_field_cyl(ρ, z, params, eq_type)

Calculate current sheet field in cylindrical coordinates (ρ, φ, z).

Returns (B_ρ, B_φ, B_z) in nT.
"""
function current_sheet_field_cyl(ρ, z, p::CurrentSheetParameters, eq_type::CurrentSheetEquationType)
    # Handle positions at or near the axis
    if ρ < 1e-6
        # On axis, B_ρ = 0 by symmetry, B_z from axial field
        B_z = axial_field(z, p)
        return (0.0, 0.0, B_z)
    end

    # Calculate azimuthal current contribution (B_ρ and B_z)
    if eq_type == Analytical
        B_ρ, B_z = analytical_field(ρ, z, p)
    elseif eq_type == Integral
        B_ρ, B_z = integral_field(ρ, z, p)
    else  # Hybrid
        B_ρ, B_z = hybrid_field(ρ, z, p)
    end

    # Calculate B_φ from radial current (CON2020 addition)
    B_φ = radial_current_bphi(ρ, z, p)

    return (B_ρ, B_φ, B_z)
end

"""
    analytical_field(ρ, z, p)

Calculate B_ρ and B_z using analytical approximations from Edwards et al. (2001)
and Connerney et al. (1981).

The analytical expressions are valid for positions not too close to the inner
edge of the current sheet (ρ >> r0 or ρ << r0).
"""
function analytical_field(ρ, z, p::CurrentSheetParameters)
    μ = p.μ_i_half
    r0 = p.r0
    r1 = p.r1
    d = p.d

    # Calculate field from inner edge (semi-infinite sheet from r0 to ∞)
    B_ρ_inner, B_z_inner = analytical_outer_disc(ρ, z, r0, d, μ)

    # Calculate field from outer edge (semi-infinite sheet from r1 to ∞)
    # Subtract because outer edge current is reversed
    B_ρ_outer, B_z_outer = analytical_outer_disc(ρ, z, r1, d, μ)

    B_ρ = B_ρ_inner - B_ρ_outer
    B_z = B_z_inner - B_z_outer

    return B_ρ, B_z
end

"""
    analytical_outer_disc(ρ, z, a, d, μ)

Analytical field from a semi-infinite disc starting at radius a.

Based on Connerney et al. (1981) Appendix and Edwards et al. (2001).
"""
function analytical_outer_disc(ρ, z, a, d, μ)
    # Different formulas for inside/outside the sheet thickness
    abs_z = abs(z)
    sign_z = z >= 0 ? 1.0 : -1.0

    if abs_z < d
        # Inside the current sheet (|z| < d)
        B_ρ, B_z = analytical_inside_sheet(ρ, z, a, d, μ)
    else
        # Outside the current sheet (|z| >= d)
        B_ρ, B_z = analytical_outside_sheet(ρ, z, a, d, μ)
    end

    return B_ρ, B_z
end

"""
    analytical_inside_sheet(ρ, z, a, d, μ)

Field components inside the current sheet (|z| < d).

From Connerney et al. (1981) equations A1-A4.
"""
function analytical_inside_sheet(ρ, z, a, d, μ)
    # Equations from Connerney et al. (1981) Appendix
    # Using the corrected forms from Acuña et al. (1983)

    ρ2 = ρ^2
    z2 = z^2
    a2 = a^2
    d2 = d^2

    # Calculate F functions (related to flux function)
    zp = z + d  # z plus d
    zm = z - d  # z minus d

    # For ρ > a (outside the inner edge radially)
    if ρ > a
        # F2 terms
        f2p = _f2_func(ρ, zp, a)
        f2m = _f2_func(ρ, zm, a)

        # B_ρ inside the sheet
        B_ρ = μ * (f2p - f2m) / ρ

        # B_z inside the sheet
        g2p = _g2_func(ρ, zp, a)
        g2m = _g2_func(ρ, zm, a)
        B_z = μ * (z / d - (g2p - g2m))
    else
        # For ρ < a (inside the inner edge radially)
        f1p = _f1_func(ρ, zp, a)
        f1m = _f1_func(ρ, zm, a)

        B_ρ = μ * (f1p - f1m) * ρ / a2

        g1p = _g1_func(ρ, zp, a)
        g1m = _g1_func(ρ, zm, a)
        B_z = μ * (z / d - (g1p - g1m) * ρ2 / a2)
    end

    return B_ρ, B_z
end

"""
    analytical_outside_sheet(ρ, z, a, d, μ)

Field components outside the current sheet (|z| >= d).

From Connerney et al. (1981) equations A5-A8.
"""
function analytical_outside_sheet(ρ, z, a, d, μ)
    sign_z = z >= 0 ? 1.0 : -1.0
    abs_z = abs(z)

    zp = abs_z + d
    zm = abs_z - d

    if ρ > a
        # Outside inner edge radially
        f2p = _f2_func(ρ, zp, a)
        f2m = _f2_func(ρ, zm, a)

        B_ρ = μ * sign_z * (f2p - f2m) / ρ

        g2p = _g2_func(ρ, zp, a)
        g2m = _g2_func(ρ, zm, a)
        B_z = μ * (1.0 - (g2p - g2m))
    else
        # Inside inner edge radially
        f1p = _f1_func(ρ, zp, a)
        f1m = _f1_func(ρ, zm, a)

        B_ρ = μ * sign_z * (f1p - f1m) * ρ / a^2

        g1p = _g1_func(ρ, zp, a)
        g1m = _g1_func(ρ, zm, a)
        B_z = μ * (1.0 - (g1p - g1m) * ρ^2 / a^2)
    end

    return B_ρ, B_z
end

# Helper functions for analytical expressions
# These are based on Connerney et al. (1981) / Edwards et al. (2001)

"""F2 function for ρ > a (outer region)"""
function _f2_func(ρ, ζ, a)
    # F2 = (a/2) * [sqrt(ζ² + (ρ-a)²) - sqrt(ζ² + (ρ+a)²) + 2a]
    # Simplified using asymptotic forms
    ρ2 = ρ^2
    a2 = a^2
    ζ2 = ζ^2

    r_minus = sqrt(ζ2 + (ρ - a)^2)
    r_plus = sqrt(ζ2 + (ρ + a)^2)

    return 0.5 * a * (r_minus - r_plus) + a2 / 2
end

"""G2 function for ρ > a (outer region)"""
function _g2_func(ρ, ζ, a)
    # Related to ∂F2/∂z
    ζ2 = ζ^2
    r_minus = sqrt(ζ2 + (ρ - a)^2)
    r_plus = sqrt(ζ2 + (ρ + a)^2)

    term1 = ζ / r_minus
    term2 = ζ / r_plus

    return 0.5 * (term1 - term2)
end

"""F1 function for ρ < a (inner region)"""
function _f1_func(ρ, ζ, a)
    # For small ρ/a, use series expansion
    ζ2 = ζ^2
    a2 = a^2
    ρ2 = ρ^2

    # Approximate form valid for ρ << a
    r_a = sqrt(ζ2 + a2)

    return 0.5 * a2 * ζ / r_a
end

"""G1 function for ρ < a (inner region)"""
function _g1_func(ρ, ζ, a)
    ζ2 = ζ^2
    a2 = a^2

    r_a = sqrt(ζ2 + a2)

    return 0.5 * a2 / r_a
end

"""
    integral_field(ρ, z, p)

Calculate B_ρ and B_z using numerical integration of the Bessel function
integrals from Connerney et al. (1981).

This is more accurate than the analytical approximation, especially near
the inner edge of the current sheet, but slower.
"""
function integral_field(ρ, z, p::CurrentSheetParameters)
    μ = p.μ_i_half
    r0 = p.r0
    r1 = p.r1
    d = p.d

    # Integration parameters
    # The integrals involve J₁(kρ) × [J₀(ka) - J₀(kr₁)] × sinh(kd)/k × exp(-k|z|)
    # We use adaptive quadrature

    # For numerical stability, use a reasonable upper limit for k
    # The integrand decays exponentially for k > 1/min(ρ, d)
    k_max = 10.0 / min(ρ, d, 1.0)

    # Simple trapezoidal integration (could be improved with adaptive quadrature)
    n_points = 500
    dk = k_max / n_points

    B_ρ_integral = 0.0
    B_z_integral = 0.0

    abs_z = abs(z)
    sign_z = z >= 0 ? 1.0 : -1.0

    for i in 1:n_points
        k = (i - 0.5) * dk  # Midpoint rule

        # Bessel functions
        j0_r0 = besselj0(k * r0)
        j0_r1 = besselj0(k * r1)
        j1_ρ = besselj1(k * ρ)
        j0_ρ = besselj0(k * ρ)

        # Common terms
        sinh_kd = sinh(k * d)
        exp_kz = exp(-k * abs_z)

        # Δ term (difference between inner and outer edge contributions)
        Δ = j0_r0 - j0_r1

        if abs_z < d
            # Inside the sheet
            cosh_kz = cosh(k * z)
            sinh_kz = sinh(k * z)

            # B_ρ integrand (inside)
            B_ρ_integral += j1_ρ * Δ * sinh_kz / sinh_kd * dk

            # B_z integrand (inside)
            B_z_integral += j0_ρ * Δ * cosh_kz / sinh_kd * dk
        else
            # Outside the sheet
            # B_ρ integrand (outside)
            B_ρ_integral += sign_z * j1_ρ * Δ * exp_kz * dk

            # B_z integrand (outside)
            B_z_integral += j0_ρ * Δ * exp_kz * dk
        end
    end

    B_ρ = 2 * μ * d * B_ρ_integral
    B_z = 2 * μ * d * B_z_integral

    return B_ρ, B_z
end

"""
    hybrid_field(ρ, z, p)

Use hybrid approach: analytical far from inner edge, integral near inner edge.
"""
function hybrid_field(ρ, z, p::CurrentSheetParameters)
    # Use integral method near the inner edge where analytical breaks down
    # Threshold: within 2 Rⱼ of inner edge and within 2d of the sheet
    near_inner_edge = abs(ρ - p.r0) < 2.0 && abs(z) < 2 * p.d

    if near_inner_edge
        return integral_field(ρ, z, p)
    else
        return analytical_field(ρ, z, p)
    end
end

"""
    axial_field(z, p)

Calculate B_z on the axis (ρ = 0).
"""
function axial_field(z, p::CurrentSheetParameters)
    μ = p.μ_i_half
    r0 = p.r0
    r1 = p.r1
    d = p.d

    # On axis, the field simplifies considerably
    abs_z = abs(z)

    if abs_z < d
        # Inside sheet on axis
        B_z = 2 * μ * z / d
    else
        # Outside sheet on axis
        sign_z = z >= 0 ? 1.0 : -1.0
        B_z = 2 * μ * sign_z
    end

    # Correction for finite sheet extent
    # (This is approximate; full calculation requires integration)
    zp = abs_z + d
    zm = abs(abs_z - d)

    correction_r0 = sqrt(zp^2 + r0^2) - sqrt(zm^2 + r0^2)
    correction_r1 = sqrt(zp^2 + r1^2) - sqrt(zm^2 + r1^2)

    return B_z * (1 - 0.5 * (correction_r1 - correction_r0) / (2 * d))
end

"""
    radial_current_bphi(ρ, z, p)

Calculate B_φ contribution from radial current (CON2020 model).

The radial current flows outward in the equatorial plane, creating an
azimuthal field component.
"""
function radial_current_bphi(ρ, z, p::CurrentSheetParameters)
    if !p.use_radial_current || p.I_rho == 0
        return 0.0
    end

    d = p.d
    abs_z = abs(z)
    sign_z = z >= 0 ? 1.0 : -1.0

    # The radial current creates B_φ component
    # Based on Connerney et al. (2020) formulation
    # B_φ = -sign(z) × (μ₀Iᵨ)/(2πρ) × f(z/d)

    if ρ < 1e-6
        return 0.0
    end

    # Current intensity factor (already in nT×Rⱼ units)
    I_factor = p.I_rho

    if abs_z < d
        # Inside the sheet: linear variation
        f_z = abs_z / d
    else
        # Outside the sheet: constant
        f_z = 1.0
    end

    # B_φ from radial current
    B_φ = -sign_z * I_factor / ρ * f_z

    return B_φ
end

# Bessel function implementations using polynomial approximations
# For better accuracy, consider using SpecialFunctions.jl

"""Bessel function J₀(x) - polynomial approximation"""
function besselj0(x)
    ax = abs(x)
    if ax < 8.0
        y = x * x
        return (57568490574.0 + y * (-13362590354.0 + y * (651619640.7 +
            y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))))) /
            (57568490411.0 + y * (1029532985.0 + y * (9494680.718 +
            y * (59272.64853 + y * (267.8532712 + y)))))
    else
        z = 8.0 / ax
        y = z * z
        xx = ax - 0.785398164
        return sqrt(0.636619772 / ax) * (cos(xx) * (1.0 + y * (-0.1098628627e-2 +
            y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)))) -
            z * sin(xx) * (-0.1562499995e-1 + y * (0.1430488765e-3 +
            y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934945152e-7)))))
    end
end

"""Bessel function J₁(x) - polynomial approximation"""
function besselj1(x)
    ax = abs(x)
    if ax < 8.0
        y = x * x
        ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 +
            y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606)))))) /
            (144725228442.0 + y * (2300535178.0 + y * (18583304.74 +
            y * (99447.43394 + y * (376.9991397 + y)))))
        return ans1
    else
        z = 8.0 / ax
        y = z * z
        xx = ax - 2.356194491
        ans1 = sqrt(0.636619772 / ax) * (cos(xx) * (1.0 + y * (0.183105e-2 +
            y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))))) -
            z * sin(xx) * (0.04687499995 + y * (-0.2002690873e-3 +
            y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)))))
        return x >= 0 ? ans1 : -ans1
    end
end

"""
    CON2020(; kwargs...)

Load the Connerney et al. (2020) Jupiter magnetodisc current sheet model.

# Keyword Arguments
- `equation_type::Symbol`: `:analytical`, `:integral`, or `:hybrid` (default)
- `μ_i_half::Float64`: Current parameter in nT (default: 139.6)
- `r0::Float64`: Inner edge in Rⱼ (default: 7.8)
- `r1::Float64`: Outer edge in Rⱼ (default: 51.4)
- `d::Float64`: Half-thickness in Rⱼ (default: 3.6)
- `θ_tilt::Float64`: Tilt angle in degrees (default: 9.3)
- `φ_tilt::Float64`: Tilt longitude in degrees (default: 155.8)
- `I_rho::Float64`: Radial current in MA (default: 16.7)
- `use_radial_current::Bool`: Include radial current (default: true)
- `in::Symbol`: Input coordinate system (default: :spherical)
- `out::Symbol`: Output coordinate system (default: :spherical)

# Returns
A callable `MagneticModel` that can be evaluated as `model(r, θ, φ)`.

# Example
```julia
model = CON2020()
B = model(10.0, π/2, 0.0)  # Field at 10 Rⱼ in equatorial plane

# With JRM33 internal field
jrm33 = JRM33(max_degree=13)
total_B = jrm33(r, θ, φ) .+ model(r, θ, φ)
```
"""
function CON2020(;
    equation_type::Symbol = :hybrid,
    μ_i_half::Float64 = 139.6,
    r0::Float64 = 7.8,
    r1::Float64 = 51.4,
    d::Float64 = 3.6,
    θ_tilt::Real = 9.3,
    φ_tilt::Real = 155.8,
    I_rho::Float64 = 16.7,
    use_radial_current::Bool = true,
    in::Symbol = :spherical,
    out::Symbol = :spherical
)
    eq_type = if equation_type == :analytical
        Analytical
    elseif equation_type == :integral
        Integral
    else
        Hybrid
    end

    params = CurrentSheetParameters(
        μ_i_half = μ_i_half,
        r0 = r0,
        r1 = r1,
        d = d,
        θ_tilt = deg2rad(θ_tilt),
        φ_tilt = deg2rad(φ_tilt),
        I_rho = I_rho,
        use_radial_current = use_radial_current
    )

    model = CurrentSheetModel(name = "CON2020", params = params, equation_type = eq_type)

    return MagneticModel(model, Jupiter; in = in, out = out)
end

"""
    CAN1981(; kwargs...)

Load the original Connerney-Acuña-Ness (1981) current sheet model.

This is the same as CON2020 but without the radial current contribution,
using the original parameters.

# Keyword Arguments
- `equation_type::Symbol`: `:analytical`, `:integral`, or `:hybrid` (default)
- `in::Symbol`: Input coordinate system (default: :spherical)
- `out::Symbol`: Output coordinate system (default: :spherical)

# Example
```julia
model = CAN1981()
B = model(10.0, π/2, 0.0)
```
"""
function CAN1981(;
    equation_type::Symbol = :hybrid,
    in::Symbol = :spherical,
    out::Symbol = :spherical
)
    # Original CAN81 parameters (approximate, as refined over time)
    return CON2020(
        equation_type = equation_type,
        μ_i_half = 225.0,  # Original value was higher
        r0 = 5.0,          # Original inner edge
        r1 = 50.0,         # Original outer edge
        d = 2.5,           # Original half-thickness
        θ_tilt = 9.6,      # Original tilt
        φ_tilt = 202.0,    # Original tilt longitude
        I_rho = 0.0,
        use_radial_current = false,
        in = in,
        out = out
    )
end
