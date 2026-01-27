"""
Core type definitions for PlanetaryMagneticFields.jl
"""

struct Planet
    name::Symbol
    radius::Float64  # Mean radius in km
end

struct PlanetFrame <: AbstractReferenceFrame end

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

evalsph(coeffs::GaussCoefficients, r, θ, φ) = evalsph(coeffs, r, θ, φ, degree(coeffs), order(coeffs))


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
    t = eltype(eps)(t)

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
- `obj`: Object (planet, body) for which the model is defined
- `degree`: Maximum degree N for evaluation
- `order`: Maximum order M for evaluation
"""
struct SphericalHarmonicModel{C, O, F} <: InternalFieldModel
    name::String
    coeffs::C
    obj::O
    degree::Int
    order::Int
    frame::F
end

function SphericalHarmonicModel(name, coeffs, obj = nothing; degree = nothing, order = nothing, frame = PlanetFrame())
    degree = @something degree PlanetaryMagneticFields.degree(coeffs)
    order = min(@something(order, PlanetaryMagneticFields.order(coeffs)), degree)
    return SphericalHarmonicModel(name, coeffs, obj, degree, order, frame)
end

getcsys(m::SphericalHarmonicModel) = (m.frame, Spherical())

# Convenience accessors
degree(model::SphericalHarmonicModel) = model.degree
order(model::SphericalHarmonicModel) = model.order
epochs(model::SphericalHarmonicModel{<:TimeVaryingGaussCoefficients}) = model.coeffs.epochs
epoch_range(model::SphericalHarmonicModel{<:TimeVaryingGaussCoefficients}) = (model.coeffs.epochs[1], model.coeffs.epochs[end])

# Check if a model is time-varying
is_time_varying(::SphericalHarmonicModel{<:GaussCoefficients}) = false
is_time_varying(::SphericalHarmonicModel{<:TimeVaryingGaussCoefficients}) = true

include("show.jl")
