# Define a GF2^16 element of a binary field
struct BinaryElem16 <: BinaryElem
    poly::BinaryPoly16
end

irreducible_poly(::Type{BinaryElem16}) = BinaryPoly16(UInt16(0b101101)) # a^16 + a^5 + a^3 + a^2 + 1, standard
# Make these macro-generated
function mod_irreducible(a::BinaryPoly32)
    (hi, lo) = split(a)

    tmp = hi + (hi >> (16 - 5)) + (hi >> (16 - 3)) + (hi >> (16 - 2))
    res = lo + tmp + (tmp << 2) + (tmp << 3) + (tmp << 5)

    return BinaryElem16(res)
end

Elem2x16 = NTuple{2, BinaryElem16}

# Hopefully this automatically SIMDs
function *(a::Elem2x16, b::Elem2x16)
    mod_irreducible.(poly.(a) * poly.(b))
end

function *(λ::BinaryElem16, b::Elem2x16)
    λp = poly(λ)
    mod_irreducible.(λp * poly.(b))
end

+(a::Elem2x16, b::Elem2x16) = a .+ b