# Define a GF2^128 element of a binary field
struct BinaryElem128 <: BinaryElem
    poly::BinaryPoly128
end
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{BinaryElem128}) = BinaryElem128(rand(rng, BinaryPoly128))
Base.convert(::Type{BinaryElem128}, v::BinaryPoly128) = BinaryElem128(v)
Base.convert(::Type{BinaryElem128}, v::UInt128) = BinaryElem128(BinaryPoly128(v))

# tmp hack 
Base.zero(::Type{BinaryElem128}) = BinaryPoly128(UInt128(0))
Base.one(::Type{BinaryElem128}) = BinaryPoly128(UInt128(1))


irreducible_poly(::Type{BinaryElem128}) = BinaryPoly128(UInt128(0b10000111)) # x^128 + x^7 + x^2 + x + 1, standard

# Make these macro-generated
function mod_irreducible(a::NTuple{2, BinaryPoly128})
    (hi, lo) = a

    tmp = hi + (hi >> 127) + (hi >> 126) + (hi >> 121)
    res = lo + tmp + (tmp << 1) + (tmp << 2) + (tmp << 7)

    return convert(BinaryElem128, res)
end
