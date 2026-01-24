"""
Core type definitions for MagneticModels.jl

This module defines the abstract type hierarchy and concrete types for representing
magnetic field models.
"""

"""
    MagneticFieldModel

Abstract base type for all magnetic field models.
"""
abstract type MagneticFieldModel end

"""
    InternalFieldModel <: MagneticFieldModel

Abstract type for internal (planetary) magnetic field models.
These models typically use spherical harmonic expansions.
"""
abstract type InternalFieldModel <: MagneticFieldModel end

"""
    ExternalFieldModel <: MagneticFieldModel

Abstract type for external (magnetospheric) magnetic field models.
These models may be empirical or physics-based.
"""
abstract type ExternalFieldModel <: MagneticFieldModel end

struct Planet
    name::Symbol
    radius::Float64
end

const Jupiter = Planet(:jupiter, 71492.0)
const Saturn = Planet(:saturn, 60268.0)
const Earth = Planet(:earth, 6371.0)

planet(s::Symbol) = if s in (:jupiter, :Jupiter)
    Jupiter
elseif s in (:saturn, :Saturn)
    Saturn
elseif s in (:earth, :Earth)
    Earth
else
    error("Unknown planet: $s")
end

"""
    GaussCoefficients

Storage for Gauss coefficients in Schmidt semi-normalized form.

# Fields
- `g::Matrix{Float64}`: Schmidt semi-normalized g coefficients (degree × order)
- `h::Matrix{Float64}`: Schmidt semi-normalized h coefficients (degree × order)
- `degree::Int`: Maximum degree N
- `order::Int`: Maximum order M (typically M ≤ N)

# Indexing
Coefficients are stored in a triangular matrix where row n (degree) contains
coefficients for orders m = 0, 1, ..., min(n, M).

For coefficient g_n^m or h_n^m:
- Use: coeffs.g[n+1, m+1] or coeffs.h[n+1, m+1]
- Note: Julia uses 1-based indexing, so degree n is in row n+1

# Schmidt Semi-Normalization
These coefficients follow the Schmidt semi-normalization convention used in
geomagnetism and planetary magnetic field modeling (IGRF, WMM, JRM, etc.).
"""
struct GaussCoefficients
    g::Matrix{Float64}
    h::Matrix{Float64}
    degree::Int
    order::Int

    function GaussCoefficients(g::Matrix{Float64}, h::Matrix{Float64}, degree::Int, order::Int)
        # Validation
        size(g) == size(h) || error("g and h matrices must have the same size")
        size(g, 1) >= degree + 1 || error("g matrix must have at least degree+1 rows")
        size(g, 2) >= order + 1 || error("g matrix must have at least order+1 columns")
        degree >= 1 || error("Degree must be at least 1")
        order >= 0 || error("Order must be non-negative")
        order <= degree || error("Order cannot exceed degree")

        return new(g, h, degree, order)
    end
end

"""
    Base.getindex(coeffs::GaussCoefficients, n::Int, m::Int)

Get the (g, h) coefficient pair for degree n and order m.

Returns a named tuple with fields `g` and `h`.

# Example
```julia
coeffs[2, 1]  # Get g_2^1 and h_2^1
```
"""
function Base.getindex(coeffs::GaussCoefficients, n::Int, m::Int)
    return (g = coeffs.g[n + 1, m + 1], h = coeffs.h[n + 1, m + 1])
end

"""
    SphericalHarmonicModel <: InternalFieldModel

A magnetic field model using spherical harmonic expansion.

# Fields
- `name::String`: Model name (e.g., "jrm09", "jrm33", "igrf14")
- `coeffs::GaussCoefficients`: Gauss coefficients

# Mathematical Formulation
The magnetic scalar potential V is expanded in spherical harmonics:

```math
V(r,θ,φ) = a \\sum_{n=1}^{N} \\sum_{m=0}^{n} \\left(\\frac{a}{r}\\right)^{n+1}
           [g_n^m \\cos(mφ) + h_n^m \\sin(mφ)] P_n^m(\\cos θ)
```

where:
- a is the planetary radius
- g_n^m, h_n^m are the Gauss coefficients (Schmidt semi-normalized)
- P_n^m are the associated Legendre polynomials (Schmidt semi-normalized)
- (r,θ,φ) are spherical coordinates (radius, colatitude, longitude)

The magnetic field **B** = -∇V is computed from derivatives of the potential.

# Coordinate Systems
- Spherical: (r, θ, φ) where θ is colatitude [0, π] and φ is longitude [0, 2π]
- Cartesian: (x, y, z) in the planet's reference frame
- All coordinates are planetocentric unless otherwise specified
"""
struct SphericalHarmonicModel <: InternalFieldModel
    name::String
    coeffs::GaussCoefficients
end

# Note: Input positions are assumed to be in planetary radii by default.
function (m::SphericalHarmonicModel)(r, θ, φ)
    return evaluate_field_spherical(m, r, θ, φ)
end

"""
    degree(model::SphericalHarmonicModel)

Get the maximum degree of the spherical harmonic expansion.
"""
degree(model::SphericalHarmonicModel) = model.coeffs.degree

"""
    order(model::SphericalHarmonicModel)

Get the maximum order of the spherical harmonic expansion.
"""
order(model::SphericalHarmonicModel) = model.coeffs.order

"""
    MagneticModel

A callable wrapper around a magnetic field model that stores default parameters.

# Usage
```julia
model = JRM33()
B = model(r, θ, φ)  # Uses default in/out from constructor
B = model(r, θ, φ; out=:spherical)  # Keyword arguments override defaults
```
"""
struct MagneticModel{M, O}
    model::M
    obj::O
    in::Symbol
    out::Symbol
end

function MagneticModel(model, obj; in = :spherical, out = :spherical)
    return MagneticModel(model, obj, in, out)
end

function Base.getproperty(m::MagneticModel, s::Symbol)
    return s in fieldnames(typeof(m)) ? getfield(m, s) : getproperty(m.model, s)
end

# Note: Input positions are assumed to be in planetary radii by default.
# For physical units (km), use Unitful.jl via the package extension.
function (m::MagneticModel)(
        r, θ, φ;
        in = nothing,
        out = nothing
    )
    # Use keyword arguments if provided, otherwise use defaults from constructor
    coords = isnothing(in) ? m.in : in
    output_coords = isnothing(out) ? m.out : out

    # Validate coordinate systems
    coords ∈ (:spherical, :cartesian) || error("in must be :spherical or :cartesian")
    output_coords ∈ (:spherical, :cartesian) || error("out must be :spherical or :cartesian")

    # Create position vector (assumed to be in planetary radii)
    pos = [Float64(r), Float64(θ), Float64(φ)]

    # Evaluate based on coordinate systems
    return if coords == :spherical && output_coords == :spherical
        m.model(r, θ, φ)
    elseif coords == :spherical && output_coords == :cartesian
        B_sph = m.model(r, θ, φ)
        spherical_field_to_cartesian(B_sph..., pos[2], pos[3])
    elseif coords == :cartesian && output_coords == :spherical
        (r_sph, θ_sph, φ_sph) = cartesian_to_spherical(pos...)
        m.model(r_sph, θ_sph, φ_sph)
    else  # cartesian -> cartesian
        evaluate_field_cartesian(m.model, pos...)
    end
end

(m::MagneticModel)(r) = m(r...)


# Convenience accessors
degree(m::MagneticModel) = degree(m.model)
order(m::MagneticModel) = order(m.model)

include("show.jl")
