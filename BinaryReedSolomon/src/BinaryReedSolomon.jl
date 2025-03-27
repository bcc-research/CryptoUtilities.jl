module BinaryReedSolomon

using BinaryFields

# Utility functions
is_pow_2(n) = 2^(round(Int, log2(n))) == n


export fft!, ifft!, compute_twiddles, compute_twiddles!
include("./binaryfft.jl")

export reed_solomon, encode, encode!
include("./reedsolomon.jl")

end # module BinaryReedSolomon
