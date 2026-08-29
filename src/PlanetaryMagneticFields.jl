module PlanetaryMagneticFields

using LinearAlgebra
using TOML
using Dates
using GeoCotrans: Cartesian3, Spherical, sph2car, car2sph, GEO, InternalFieldModel, CompositeFieldModel, AbstractReferenceFrame, IGRF, LinearInterp
import GeoCotrans: getcsys, evalsph, evalmodel, trace

include("types.jl")
include("planets.jl")
include("data.jl")
include("coefficients.jl")
include("spherical_harmonics.jl")
include("api.jl")

export load_model, available_models, model_info, fieldmap
export trace

# Convenience model accessors
export JRM09,
    JRM33,
    IGRF

# Utility functions that users might want
export degree, order
export getcsys, Cartesian3, Spherical

export plot_fieldmap, plot_fieldmap!, plot_models

end
