module BinaryFields

using SIMD, Random
import Base: +, *, zero, transpose, adjoint, inv, <<, >>

export BinaryPoly8, BinaryPoly16, BinaryPoly32, BinaryPoly64, BinaryPoly128
export BinaryElem, BinaryElem8, BinaryElem16, BinaryElem32, BinaryElem64, BinaryElem128

include("./utilities.jl")

include("./binarypoly.jl")
include("./binaryfield.jl")

end # module BinaryFields
