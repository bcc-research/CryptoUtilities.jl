# Define a GF2^32 element of a binary field
struct BinaryElem32 <: BinaryElem
    poly::BinaryPoly32
end

irreducible_poly(::Type{BinaryElem32}) = BinaryPoly32(UInt32(0b11001 | 1 << 7 | 1 << 9 | 1 << 15)) # a^32 + a^15 + a^9 + a^7 + x^4 + x^3 + 1, standard

# XXX: Make these macro-generated
function mod_irreducible(a::BinaryPoly64)
    (hi, lo) = split(a)

    tmp = hi + (hi >> (32 - 15)) + (hi >> (32 - 9)) + (hi >> (32 - 7)) + (hi >> (32 - 4)) + (hi >> (32 - 3))
    res = lo + tmp + (tmp << 3) + (tmp << 4) + (tmp << 7) + (tmp << 9) + (tmp << 15)

    return BinaryElem32(res)
end