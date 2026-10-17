module PlanetaryMagneticFields

using LinearAlgebra
using TOML
using Dates
using GeoCotrans: Cartesian3, Spherical, sph2car, AbstractReferenceFrame, IGRF
using GeoCotrans: GaussCoefficients
import GeoCotrans: SphericalHarmonicModel, getcsys, trace

include("types.jl")
include("planets.jl")
include("data.jl")
include("coefficients.jl")
include("spherical_harmonics.jl")
include("api.jl")

export load_model, available_models, model_info, fieldmap
export trace

export JRM09, JRM33, IGRF
export degree, order
export getcsys, Cartesian3, Spherical

export plot_fieldmap, plot_fieldmap!, plot_models

end
