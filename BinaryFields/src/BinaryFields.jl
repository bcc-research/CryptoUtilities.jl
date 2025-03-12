module BinaryFields

using SIMD, Random
import Base: +, *, zero, one, transpose, adjoint, inv, <<, >>

export BinaryPoly16, BinaryPoly64, BinaryPoly128
export BinaryElem, BinaryElem16, BinaryElem128

include("./utilities.jl")

include("./binarypoly.jl")
include("./binaryfield.jl")

include("./be128.jl")
include("./be16.jl")


end # module BinaryFields
