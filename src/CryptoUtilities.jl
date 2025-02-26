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

# Other crypto stuff
include("./merkletree.jl")
include("./merkle_proofs.jl")


end # module CelestiaProve
