# Define a GF2^16 element of a binary field
struct BinaryElem16 <: BinaryElem
    poly::BinaryPoly16
end
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{BinaryElem16}) = BinaryElem16(rand(rng, BinaryPoly16))
Base.convert(::Type{BinaryElem16}, v::BinaryPoly16) = BinaryElem16(v)
Base.convert(::Type{BinaryElem16}, v::UInt16) = BinaryElem16(BinaryPoly16(v))

# tmp hack 
Base.one(::Type{BinaryElem16}) = BinaryPoly16(UInt16(1))


irreducible_poly(::Type{BinaryElem16}) = BinaryPoly16(UInt16(0b101101)) # a^16 + a^5 + a^3 + a^2 + 1, standard
# Make these macro-generated
function mod_irreducible(a::NTuple{2, BinaryPoly16})
    (hi, lo) = a

    tmp = hi + (hi >> (16 - 5)) + (hi >> (16 - 3)) + (hi >> (16 - 2))
    res = lo + tmp + (tmp << 2) + (tmp << 3) + (tmp << 5)

    return convert(BinaryElem16, res)
end
