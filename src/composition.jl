"""
    CompositeModel

Combine multiple `MagneticModel`s into a single model by summing their fields.
"""
struct CompositeModel{M<:Tuple} <: MagneticFieldModel
    models::M
end

function CompositeModel(models::Vararg{MagneticModel})
    isempty(models) && error("At least one model is required to build a CompositeModel")
    planets = getfield.(models, :obj)
    length(unique(planets)) == 1 || error("All models in a CompositeModel must reference the same planet")
    return CompositeModel(tuple(models...))
end

(m::CompositeModel)(r, θ, φ) = evaluate_field_spherical(m, r, θ, φ)

function evaluate_field_spherical(m::CompositeModel, r, θ, φ)
    return mapreduce(model -> model(r, θ, φ), +, m.models)
end

"""
    compose_models(models...; in=:spherical, out=:spherical)

Create a `MagneticModel` that represents the sum of the provided `models`.
All models must reference the same planet.
"""
function compose_models(models::Vararg{MagneticModel}; in = :spherical, out = :spherical)
    comp = CompositeModel(models...)
    return MagneticModel(comp, models[1].obj; in = in, out = out)
end

compose_models(models::Tuple; kw...) = compose_models(models...; kw...)
