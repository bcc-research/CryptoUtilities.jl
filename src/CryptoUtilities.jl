module CryptoUtilities

include("./utilities.jl")

using LinearAlgebra

# Internal packages
using BinaryFields

# Binary FFT stuff
include("./binaryfft.jl")

export reed_solomon
include("./reedsolomon.jl")

include("./experiment.jl")

# Other crypto stuff
# Commenting out for now to be able to run tests
# include("./merkletree.jl")
# include("./merkle_proofs.jl")


end # module CelestiaProve
