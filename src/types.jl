"""
Core type definitions for PlanetaryMagneticFields.jl
"""

"""
    MagneticFieldModel

Abstract base type for all magnetic field models.
"""
abstract type MagneticFieldModel end

"""
    InternalFieldModel <: MagneticFieldModel

Abstract type for internal (planetary) magnetic field models.
"""
abstract type InternalFieldModel <: MagneticFieldModel end

"""
    ExternalFieldModel <: MagneticFieldModel

Abstract type for external (magnetospheric) magnetic field models.
"""
abstract type ExternalFieldModel <: MagneticFieldModel end

struct Planet
    name::Symbol
    radius::Float64  # Mean radius in km
end

"""
    GaussCoefficients

Storage for Gauss coefficients in Schmidt semi-normalized form.

# Fields
- `g`: Schmidt semi-normalized g coefficients (degree × order)
- `h`: Schmidt semi-normalized h coefficients (degree × order)

The degree and order are inferred from matrix dimensions: `degree = size(g,1)-1`, `order = size(g,2)-1`.
"""
struct GaussCoefficients{A}
    g::A
    h::A

    function GaussCoefficients(g::A, h::A; check = false) where {A <: AbstractMatrix}
        check && begin
            size(g) == size(h) || error("g and h matrices must have the same size")
            size(g, 1) >= 2 || error("g matrix must have at least 2 rows (degree >= 1)")
            size(g, 2) >= 1 || error("g matrix must have at least 1 column (order >= 0)")
            size(g, 2) <= size(g, 1) || error("Order cannot exceed degree")
        end
        return new{A}(g, h)
    end
end

# Accessor functions for degree and order (inferred from matrix size)
degree(coeffs::GaussCoefficients) = size(coeffs.g, 1) - 1
order(coeffs::GaussCoefficients) = size(coeffs.g, 2) - 1

# Get (g, h) coefficient pair for degree n and order m
function Base.getindex(coeffs::GaussCoefficients, n::Int, m::Int)
    return (g = coeffs.g[n + 1, m + 1], h = coeffs.h[n + 1, m + 1])
end

"""
    TimeVaryingGaussCoefficients

Storage for time-varying Gauss coefficients with multiple epochs.
Coefficients are linearly interpolated between epochs.

# Fields
- `epochs`: Sorted vector of epoch years (e.g., [1900.0, 1905.0, ..., 2025.0])
- `coefficients`: Vector of GaussCoefficients, one per epoch
"""
struct TimeVaryingGaussCoefficients{T, C}
    epochs::Vector{T}
    coefficients::Vector{C}

    function TimeVaryingGaussCoefficients(
            epochs::Vector{T},
            coefficients::Vector{C};
            check = true
        ) where {T, C}
        check && begin
            length(epochs) == length(coefficients) || error("Number of epochs must match number of coefficient sets")
            length(epochs) >= 1 || error("At least one epoch is required")
            issorted(epochs) || error("Epochs must be sorted in ascending order")
            deg = degree(coefficients[1])
            ord = order(coefficients[1])
            for c in coefficients[2:end]
                degree(c) == deg || error("All coefficient sets must have the same degree")
                order(c) == ord || error("All coefficient sets must have the same order")
            end
        end
        return new{T, C}(epochs, coefficients)
    end
end

# Accessor functions for TimeVaryingGaussCoefficients
degree(tvc::TimeVaryingGaussCoefficients) = degree(tvc.coefficients[1])
order(tvc::TimeVaryingGaussCoefficients) = order(tvc.coefficients[1])

function (tvc::TimeVaryingGaussCoefficients)(t)
    eps = tvc.epochs
    coeffs = tvc.coefficients

    # Handle edge cases
    if t <= eps[1]
        return coeffs[1]
    elseif t >= eps[end]
        return coeffs[end]
    end

    idx = searchsortedlast(eps, t)
    t0, t1 = eps[idx], eps[idx + 1]
    c0, c1 = coeffs[idx], coeffs[idx + 1]
    α = (t - t0) / (t1 - t0)

    g_interp = LazyArray(@~ (1 - α) * c0.g + α * c1.g)
    h_interp = LazyArray(@~ (1 - α) * c0.h + α * c1.h)

    return GaussCoefficients(g_interp, h_interp)
end

"""
    SphericalHarmonicModel{C} <: InternalFieldModel

A magnetic field model using spherical harmonic expansion.

# Fields
- `name`: Model name
- `coeffs`: Gauss coefficients (GaussCoefficients or TimeVaryingGaussCoefficients)
- `degree`: Maximum degree N for evaluation
- `order`: Maximum order M for evaluation
"""
struct SphericalHarmonicModel{C} <: InternalFieldModel
    name::String
    coeffs::C
    degree::Int
    order::Int
end

function SphericalHarmonicModel(name, coeffs; degree = nothing, order = nothing)
    degree = @something degree PlanetaryMagneticFields.degree(coeffs)
    order = min(@something(order, PlanetaryMagneticFields.order(coeffs)), degree)
    return SphericalHarmonicModel(name, coeffs, degree, order)
end

# Static evaluation
function (m::SphericalHarmonicModel{<:GaussCoefficients})(r, θ, φ, _ = nothing)
    return evalsph(m.coeffs, r, θ, φ; max_degree = m.degree, max_order = m.order)
end

# Time-varying evaluation
function (m::SphericalHarmonicModel{<:TimeVaryingGaussCoefficients})(r, θ, φ, t)
    return evalsph(m.coeffs(t), r, θ, φ; max_degree = m.degree, max_order = m.order)
end

# Convenience accessors
degree(model::SphericalHarmonicModel) = model.degree
order(model::SphericalHarmonicModel) = model.order
epochs(model::SphericalHarmonicModel{<:TimeVaryingGaussCoefficients}) = model.coeffs.epochs
epoch_range(model::SphericalHarmonicModel{<:TimeVaryingGaussCoefficients}) = (model.coeffs.epochs[1], model.coeffs.epochs[end])

# Check if a model is time-varying
is_time_varying(::SphericalHarmonicModel{<:GaussCoefficients}) = false
is_time_varying(::SphericalHarmonicModel{<:TimeVaryingGaussCoefficients}) = true

"""
    MagneticModel

A callable wrapper around a magnetic field model that stores default coordinate parameters.

# Usage
```julia
# Static model
model = JRM33()
B = model(r, θ, φ)
B = model(r, θ, φ; out=:cartesian)

# Time-varying model
model = IGRF()
B = model(r, θ, φ, Date(2020))
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

# Static evaluation (3 positional arguments)
function (m::MagneticModel)(r, θ, φ, t = nothing; in = nothing, out = nothing)
    coords = isnothing(in) ? m.in : in
    output_coords = isnothing(out) ? m.out : out

    coords ∈ (:spherical, :cartesian) || error("in must be :spherical or :cartesian")
    output_coords ∈ (:spherical, :cartesian) || error("out must be :spherical or :cartesian")

    return if coords == :spherical && output_coords == :spherical
        m.model(r, θ, φ, t)
    elseif coords == :spherical && output_coords == :cartesian
        B_sph = m.model(r, θ, φ, t)
        spherical_field_to_cartesian(B_sph..., θ, φ)
    elseif coords == :cartesian && output_coords == :spherical
        r_sph, θ_sph, φ_sph = cartesian_to_spherical(r, θ, φ)
        m.model(r_sph, θ_sph, φ_sph, t)
    else  # cartesian -> cartesian
        r_sph, θ_sph, φ_sph = cartesian_to_spherical(r, θ, φ)
        B_sph = m.model(r_sph, θ_sph, φ_sph, t)
        spherical_field_to_cartesian(B_sph..., θ_sph, φ_sph)
    end
end

(m::MagneticModel)(r) = m(r...)

# Convenience accessors for MagneticModel
for f in (:degree, :order, :is_time_varying, :epochs, :epoch_range)
    @eval $f(m::MagneticModel) = $f(m.model)
end

include("show.jl")
