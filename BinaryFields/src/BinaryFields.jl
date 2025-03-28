module BinaryFields

using SIMD, Random
import Base: +, *, zero, one, transpose, adjoint, inv, <<, >>, convert, length, reinterpret, ^, convert, promote_rule

export BinaryPoly16, BinaryPoly32, BinaryPoly64, BinaryPoly128
export BinaryElem, BinaryElem16, BinaryElem32, BinaryElem128
export binary_val, bitsize

include("./binarypoly.jl")
include("./binaryfield.jl")

# Necessary for preventing frankenallocations coming from complicated types
include("./warmup.jl")

end # module BinaryFields
