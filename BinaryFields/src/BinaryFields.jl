module BinaryFields

using SIMD, Random
import Base: +, *, zero, one, transpose, adjoint, inv, <<, >>, convert

export BinaryPoly16, BinaryPoly64, BinaryPoly128
export BinaryElem, BinaryElem16, BinaryElem128
export binary_val

include("./utilities.jl")

include("./binarypoly.jl")
include("./binaryfield.jl")

include("./be128.jl")
include("./be16.jl")

# Necessary for prevernting frankenallocations coming from complicated types
include("./warmup.jl")

end # module BinaryFields
