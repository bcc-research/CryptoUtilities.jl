# Define a GF2^128 element of a binary field
struct BinaryElem128 <: BinaryElem
    poly::BinaryPoly128
end

irreducible_poly(::Type{BinaryElem128}) = BinaryPoly128(UInt128(0b10000111)) # x^128 + x^7 + x^2 + x + 1, standard

# Make these macro-generated
function mod_irreducible(a::BinaryPoly256)
    (hi, lo) = split(a)

    tmp = hi + (hi >> 127) + (hi >> 126) + (hi >> 121)
    res = lo + tmp + (tmp << 1) + (tmp << 2) + (tmp << 7)

    return BinaryElem128(res)
end
