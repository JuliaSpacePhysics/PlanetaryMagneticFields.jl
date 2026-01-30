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

Enum for selecting the equation type used in current sheet calculations.
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

"""
    CurrentSheetModel <: InternalFieldModel

Connerney et al. (1981, 2020) Jupiter magnetodisc current sheet model.

This model calculates the magnetic field contribution from a washer-shaped
azimuthal current sheet centered on Jupiter's magnetic equator.
"""
struct CurrentSheetModel <: InternalFieldModel
    name::String
    params::CurrentSheetParameters
    equation_type::CurrentSheetEquationType
    obj::Planet
end

function CurrentSheetModel(;
    name::String = "CON2020",
    params::CurrentSheetParameters = CurrentSheetParameters(),
    equation_type::CurrentSheetEquationType = Hybrid
)
    return CurrentSheetModel(name, params, equation_type, Jupiter)
end

# GeoCotrans interface
getcsys(::CurrentSheetModel) = (PlanetFrame(), Spherical())

"""
    evalmodel(model::CurrentSheetModel, r, θ, φ, t)

Evaluate the current sheet magnetic field at the given position.
Returns [Br, Bθ, Bφ] in nT (SIII spherical coordinates).
"""
function evalmodel(model::CurrentSheetModel, r, θ, φ, t = nothing)
    p = model.params

    # Convert spherical SIII to Cartesian SIII
    sinθ, cosθ = sincos(θ)
    sinφ, cosφ = sincos(φ)
    x = r * sinθ * cosφ
    y = r * sinθ * sinφ
    z = r * cosθ

    # Transform to current sheet (CS) coordinates (tilted disc frame)
    x_cs, y_cs, z_cs = siii_to_cs(x, y, z, p.θ_tilt, p.φ_tilt)

    # Calculate cylindrical coordinates in CS frame
    ρ = sqrt(x_cs^2 + y_cs^2)
    φ_cs = atan(y_cs, x_cs)

    # Calculate field in cylindrical CS coordinates
    B_ρ, B_φ_cs, B_z = current_sheet_field_cyl(ρ, z_cs, p, model.equation_type)

    # Convert cylindrical to Cartesian in CS frame
    sinφ_cs, cosφ_cs = sincos(φ_cs)
    Bx_cs = B_ρ * cosφ_cs - B_φ_cs * sinφ_cs
    By_cs = B_ρ * sinφ_cs + B_φ_cs * cosφ_cs
    Bz_cs = B_z

    # Transform back to SIII Cartesian
    Bx, By, Bz = cs_to_siii(Bx_cs, By_cs, Bz_cs, p.θ_tilt, p.φ_tilt)

    # Convert Cartesian to spherical field components
    Br = Bx * sinθ * cosφ + By * sinθ * sinφ + Bz * cosθ
    Bθ = Bx * cosθ * cosφ + By * cosθ * sinφ - Bz * sinθ
    Bφ = -Bx * sinφ + By * cosφ

    return SVector{3,Float64}(Br, Bθ, Bφ)
end

# Coordinate transformations between SIII and tilted current sheet frame

function siii_to_cs(x, y, z, θ_tilt, φ_tilt)
    sin_tilt, cos_tilt = sincos(θ_tilt)
    sin_φtilt, cos_φtilt = sincos(φ_tilt)

    # Rotate around z by -φ_tilt
    x1 = x * cos_φtilt + y * sin_φtilt
    y1 = -x * sin_φtilt + y * cos_φtilt

    # Rotate around y by -θ_tilt
    x_cs = x1 * cos_tilt + z * sin_tilt
    z_cs = -x1 * sin_tilt + z * cos_tilt

    return x_cs, y1, z_cs
end

function cs_to_siii(Bx_cs, By_cs, Bz_cs, θ_tilt, φ_tilt)
    sin_tilt, cos_tilt = sincos(θ_tilt)
    sin_φtilt, cos_φtilt = sincos(φ_tilt)

    # Rotate around y by +θ_tilt
    x1 = Bx_cs * cos_tilt - Bz_cs * sin_tilt
    z1 = Bx_cs * sin_tilt + Bz_cs * cos_tilt

    # Rotate around z by +φ_tilt
    Bx = x1 * cos_φtilt - By_cs * sin_φtilt
    By = x1 * sin_φtilt + By_cs * cos_φtilt

    return Bx, By, z1
end

# Field calculation in cylindrical coordinates

function current_sheet_field_cyl(ρ, z, p::CurrentSheetParameters, eq_type::CurrentSheetEquationType)
    if ρ < 1e-6
        B_z = axial_field(z, p)
        return (0.0, 0.0, B_z)
    end

    B_ρ, B_z = if eq_type == Analytical
        analytical_field(ρ, z, p)
    elseif eq_type == Integral
        integral_field(ρ, z, p)
    else  # Hybrid
        hybrid_field(ρ, z, p)
    end

    B_φ = radial_current_bphi(ρ, z, p)

    return (B_ρ, B_φ, B_z)
end

# Analytical field calculation (Edwards et al., 2001)

function analytical_field(ρ, z, p::CurrentSheetParameters)
    B_ρ_inner, B_z_inner = analytical_outer_disc(ρ, z, p.r0, p.d, p.μ_i_half)
    B_ρ_outer, B_z_outer = analytical_outer_disc(ρ, z, p.r1, p.d, p.μ_i_half)
    return B_ρ_inner - B_ρ_outer, B_z_inner - B_z_outer
end

function analytical_outer_disc(ρ, z, a, d, μ)
    abs_z = abs(z)
    abs_z < d ? analytical_inside_sheet(ρ, z, a, d, μ) : analytical_outside_sheet(ρ, z, a, d, μ)
end

function analytical_inside_sheet(ρ, z, a, d, μ)
    zp, zm = z + d, z - d

    if ρ > a
        f2p, f2m = _f2_func(ρ, zp, a), _f2_func(ρ, zm, a)
        g2p, g2m = _g2_func(ρ, zp, a), _g2_func(ρ, zm, a)
        B_ρ = μ * (f2p - f2m) / ρ
        B_z = μ * (z / d - (g2p - g2m))
    else
        f1p, f1m = _f1_func(ρ, zp, a), _f1_func(ρ, zm, a)
        g1p, g1m = _g1_func(ρ, zp, a), _g1_func(ρ, zm, a)
        B_ρ = μ * (f1p - f1m) * ρ / a^2
        B_z = μ * (z / d - (g1p - g1m) * ρ^2 / a^2)
    end

    return B_ρ, B_z
end

function analytical_outside_sheet(ρ, z, a, d, μ)
    sign_z = z >= 0 ? 1.0 : -1.0
    abs_z = abs(z)
    zp, zm = abs_z + d, abs_z - d

    if ρ > a
        f2p, f2m = _f2_func(ρ, zp, a), _f2_func(ρ, zm, a)
        g2p, g2m = _g2_func(ρ, zp, a), _g2_func(ρ, zm, a)
        B_ρ = μ * sign_z * (f2p - f2m) / ρ
        B_z = μ * (1.0 - (g2p - g2m))
    else
        f1p, f1m = _f1_func(ρ, zp, a), _f1_func(ρ, zm, a)
        g1p, g1m = _g1_func(ρ, zp, a), _g1_func(ρ, zm, a)
        B_ρ = μ * sign_z * (f1p - f1m) * ρ / a^2
        B_z = μ * (1.0 - (g1p - g1m) * ρ^2 / a^2)
    end

    return B_ρ, B_z
end

# Helper functions for analytical expressions
_f2_func(ρ, ζ, a) = 0.5a * (sqrt(ζ^2 + (ρ - a)^2) - sqrt(ζ^2 + (ρ + a)^2)) + 0.5a^2
_g2_func(ρ, ζ, a) = 0.5 * (ζ / sqrt(ζ^2 + (ρ - a)^2) - ζ / sqrt(ζ^2 + (ρ + a)^2))
_f1_func(ρ, ζ, a) = 0.5a^2 * ζ / sqrt(ζ^2 + a^2)
_g1_func(ρ, ζ, a) = 0.5a^2 / sqrt(ζ^2 + a^2)

# Numerical integration (Connerney et al., 1981)

function integral_field(ρ, z, p::CurrentSheetParameters)
    μ, r0, r1, d = p.μ_i_half, p.r0, p.r1, p.d
    k_max = 10.0 / min(ρ, d, 1.0)
    n_points = 500
    dk = k_max / n_points

    B_ρ_int, B_z_int = 0.0, 0.0
    abs_z = abs(z)
    sign_z = z >= 0 ? 1.0 : -1.0

    for i in 1:n_points
        k = (i - 0.5) * dk
        j0_r0, j0_r1 = besselj0(k * r0), besselj0(k * r1)
        j1_ρ, j0_ρ = besselj1(k * ρ), besselj0(k * ρ)
        Δ = j0_r0 - j0_r1

        if abs_z < d
            cosh_kz, sinh_kz = cosh(k * z), sinh(k * z)
            sinh_kd = sinh(k * d)
            B_ρ_int += j1_ρ * Δ * sinh_kz / sinh_kd * dk
            B_z_int += j0_ρ * Δ * cosh_kz / sinh_kd * dk
        else
            exp_kz = exp(-k * abs_z)
            B_ρ_int += sign_z * j1_ρ * Δ * exp_kz * dk
            B_z_int += j0_ρ * Δ * exp_kz * dk
        end
    end

    return 2μ * d * B_ρ_int, 2μ * d * B_z_int
end

function hybrid_field(ρ, z, p::CurrentSheetParameters)
    near_inner_edge = abs(ρ - p.r0) < 2.0 && abs(z) < 2p.d
    near_inner_edge ? integral_field(ρ, z, p) : analytical_field(ρ, z, p)
end

function axial_field(z, p::CurrentSheetParameters)
    μ, d = p.μ_i_half, p.d
    abs_z = abs(z)
    sign_z = z >= 0 ? 1.0 : -1.0
    abs_z < d ? 2μ * z / d : 2μ * sign_z
end

function radial_current_bphi(ρ, z, p::CurrentSheetParameters)
    !p.use_radial_current && return 0.0
    ρ < 1e-6 && return 0.0

    sign_z = z >= 0 ? 1.0 : -1.0
    f_z = abs(z) < p.d ? abs(z) / p.d : 1.0
    return -sign_z * p.I_rho / ρ * f_z
end

# Bessel functions (polynomial approximations)

function besselj0(x)
    ax = abs(x)
    if ax < 8.0
        y = x * x
        (57568490574.0 + y * (-13362590354.0 + y * (651619640.7 +
            y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))))) /
            (57568490411.0 + y * (1029532985.0 + y * (9494680.718 +
            y * (59272.64853 + y * (267.8532712 + y)))))
    else
        z = 8.0 / ax
        y = z * z
        xx = ax - 0.785398164
        sqrt(0.636619772 / ax) * (cos(xx) * (1.0 + y * (-0.1098628627e-2 +
            y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)))) -
            z * sin(xx) * (-0.1562499995e-1 + y * (0.1430488765e-3 +
            y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934945152e-7)))))
    end
end

function besselj1(x)
    ax = abs(x)
    if ax < 8.0
        y = x * x
        x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 +
            y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606)))))) /
            (144725228442.0 + y * (2300535178.0 + y * (18583304.74 +
            y * (99447.43394 + y * (376.9991397 + y)))))
    else
        z = 8.0 / ax
        y = z * z
        xx = ax - 2.356194491
        ans = sqrt(0.636619772 / ax) * (cos(xx) * (1.0 + y * (0.183105e-2 +
            y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))))) -
            z * sin(xx) * (0.04687499995 + y * (-0.2002690873e-3 +
            y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)))))
        x >= 0 ? ans : -ans
    end
end

# Public API functions

"""
    CON2020(; kwargs...)

Create the Connerney et al. (2020) Jupiter magnetodisc current sheet model.

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

# Example
```julia
model = CON2020()
B = model(10.0, π/2, 0.0)  # Field at 10 Rⱼ in equatorial plane

# Combine with JRM33 internal field
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
    use_radial_current::Bool = true
)
    eq_type = equation_type == :analytical ? Analytical :
              equation_type == :integral ? Integral : Hybrid

    params = CurrentSheetParameters(;
        μ_i_half, r0, r1, d,
        θ_tilt = deg2rad(θ_tilt),
        φ_tilt = deg2rad(φ_tilt),
        I_rho, use_radial_current
    )

    return CurrentSheetModel("CON2020", params, eq_type, Jupiter)
end

"""
    CAN1981(; equation_type=:hybrid)

Create the original Connerney-Acuña-Ness (1981) current sheet model
(without radial current contribution).
"""
function CAN1981(; equation_type::Symbol = :hybrid)
    CON2020(;
        equation_type,
        μ_i_half = 225.0,
        r0 = 5.0,
        r1 = 50.0,
        d = 2.5,
        θ_tilt = 9.6,
        φ_tilt = 202.0,
        I_rho = 0.0,
        use_radial_current = false
    )
end
