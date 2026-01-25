module MagneticModels

using StaticArrays
using LinearAlgebra

# Include core functionality
include("types.jl")
include("data.jl")
include("coefficients.jl")
include("cotrans.jl")
include("spherical_harmonics.jl")

# Include planet-specific models
include("models/Jupiter.jl")

# Include public API
include("api.jl")

# Export API
# export MagneticFieldModel, InternalFieldModel, ExternalFieldModel, SphericalHarmonicModel

# Main user-facing functions
export load_model, available_models, MagneticModel

# Convenience model accessors
export JRM09,
    JRM33

# Utility functions that users might want
export degree, order

end # module MagneticModels
