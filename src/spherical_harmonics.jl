"""
Spherical harmonic field evaluation for PlanetaryMagneticFields.jl

This module implements Schmidt semi-normalized associated Legendre polynomials
and spherical harmonic magnetic field evaluation.
"""

using LinearAlgebra
using Bumper
using AxisKeys
using SatelliteToolboxLegendre: legendre!

"""
    evaluate_field_spherical(coeffs, r, θ, φ)

Evaluate the magnetic field at a point in spherical coordinates.

# Arguments
- `coeffs`: The magnetic field coefficients
- `r`: Radial distance in planetary radii (dimensionless)
- `θ`: Colatitude in radians [0, π]
- `φ`: East longitude in radians [0, 2π]

# Returns
- `SVector{3,Float64}`: Magnetic field vector [B_r, B_θ, B_φ] in nanoTesla

# Coordinate System
- Spherical coordinates (r, θ, φ) where:
  - r: radial distance from planet center
  - θ: colatitude (0 at north pole, π at south pole)
  - φ: east longitude (0 at prime meridian, increases eastward)

# Mathematical Formulation
The magnetic field is derived from the scalar potential:
```math
B_r = -\\frac{∂V}{∂r}, \\quad
B_θ = -\\frac{1}{r}\\frac{∂V}{∂θ}, \\quad
B_φ = -\\frac{1}{r\\sin θ}\\frac{∂V}{∂φ}
```
"""
function evaluate_field_spherical(coeffs, r, θ, φ; max_degree = coeffs.degree)
    @assert max_degree <= coeffs.degree
    r > 0 || error("Radius must be positive")
    0 <= θ <= π || error("Colatitude θ must be in [0, π]")
    T = promote_type(eltype(r), eltype(θ), eltype(φ))
    return @no_escape begin
        P = @alloc(T, max_degree + 1, max_degree + 1)
        sincos_mφs = @alloc(Tuple{T, T}, max_degree + 1)
        evaluate_field_spherical!(P, sincos_mφs, coeffs.g, coeffs.h, promote(r, θ, φ)..., max_degree, coeffs.order)
    end
end

@inline function evaluate_field_spherical!(P, sincos_mφs, G, H, r::T, θ::T, φ::T, max_degree, max_order) where {T}
    sinθ, cosθ = sincos(θ)
    for m in eachindex(sincos_mφs)
        sincos_mφs[m] = sincos((m - 1) * φ)
    end
    legendre!(Val(:schmidt), P, θ, max_degree)

    Br, Bθ, Bφ = 0.0, 0.0, 0.0

    # r is in planetary radii (dimensionless), so (1/r)^k gives the radial dependence
    # For internal field: B components use (1/r)^(n+2)
    ratio = 1.0 / r
    pow = ratio * ratio * ratio   # (1/r)^(n+2) for magnetic field
    flag = abs(sinθ) > 1.0e-10
    sinθ = flag ? sinθ : 1.0e-10
    # Sum over degrees and orders
    @inbounds for l in 1:max_degree
        for m in 0:min(l, max_order)
            g = G[l + 1, m + 1]
            h = H[l + 1, m + 1]
            sin_mφ, cos_mφ = sincos_mφs[m + 1]
            Pₗₘ = P[l + 1, m + 1]

            # Spherical harmonic term
            Y = g * cos_mφ + h * sin_mφ
            # Derivatives
            dY_dφ = m * (-g * sin_mφ + h * cos_mφ)

            # Contribution to field components
            # B_r = -∂V/∂r: factor of (n+1) from derivative of (a/r)^(n+1)
            Br += (l + 1) * pow * Y * Pₗₘ
            # B_θ = -(1/r)∂V/∂θ

            dPₗₘ = if m == l
                l * cosθ * Pₗₘ / sinθ
            else
                Pₗ₋₁ₘ = P[l, m + 1]
                s = sqrt((l + m) * (l - m))
                (l * cosθ * Pₗₘ - s * Pₗ₋₁ₘ) / sinθ
            end
            Bθ -= pow * Y * dPₗₘ
            # B_φ = -(1/(r sin θ))∂V/∂φ
            flag && (Bφ -= pow * dY_dφ * Pₗₘ)
        end
        pow *= ratio
    end
    return SVector{3, T}(Br, Bθ, flag ? Bφ / sinθ : Bφ)
end

function evaluate_field_cartesian(model, x, y, z)
    (r, θ, φ) = cartesian_to_spherical(x, y, z)
    B_sph = evaluate_field_spherical(model, r, θ, φ)
    return spherical_field_to_cartesian(B_sph..., θ, φ)
end

_field_func(idx::Int) = idx <= 3 ? x -> getindex(x, idx) : norm
_field_func(::Nothing) = identity
_field_func(f) = f

function fieldmap(model, r, nlat, nlon; idx = identity)
    func = _field_func(idx)
    # Create latitude and longitude grids
    # Latitude: -90 to 90 degrees (convert to colatitude for evaluation)
    # Longitude: -180 to 180 degrees
    lats = range(-90, 90, length = nlat)
    lons = range(-180, 180, length = nlon)

    # Preallocate field array (longitude × latitude for GeoMakie surface!)
    field = zeros(Float64, nlon, nlat)

    θs = @. deg2rad(90 - lats)'
    φs = @. deg2rad(mod(lons, 360))
    field = func.(model.(r, θs, φs; in = :spherical, out = :spherical))

    return KeyedArray(field; lon = lons, lat = lats)
end

function fieldmap(model, r = 1.0; nlat = 180, nlon = 360, kw...)
    return fieldmap(model, r, nlat, nlon; kw...)
end
