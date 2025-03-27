module BinaryReedSolomon

using BinaryFields

is_pow_2(n) = 2^(round(Int, log2(n))) == n

include("./binaryfft.jl")
include("./reedsolomon.jl")

end # module BinaryReedSolomon
