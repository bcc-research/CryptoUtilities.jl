module CryptoUtilities

include("./utilities.jl")

import Base: *, +, inv
import Random
using SIMD
using LinearAlgebra

# Finite field stuff
include("./binarypoly.jl")
include("./binaryfield.jl")

# Binary FFT stuff
include("./binaryfft.jl")
include("./reedsolomon.jl")

include("./experiment.jl")

# Other crypto stuff
# Commenting out for now to be able to run tests
# include("./merkletree.jl")
# include("./merkle_proofs.jl")


end # module CelestiaProve
