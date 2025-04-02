# Using SIMD for fast reinterpret. Julia's base `reinterpret` is insanely slow
# in this particular case likely due to (in this case unnecessary) type safety
# checks

abstract type BinaryElem <: Number end

poly(x::T) where T <: BinaryElem = x.poly
poly_type(::Type{T}) where T <: BinaryElem = fieldtypes(T)[1]
binary_val(x::T) where T <: BinaryElem = binary_val(poly(x))
bitsize(::Type{T}) where T <: BinaryElem = sizeof(T) * 8

Base.zero(::T) where T <: BinaryElem = T(0)
Base.zero(::Type{T}) where T <: BinaryElem = T(0)
Base.one(::T) where T <: BinaryElem = T(1)
Base.one(::Type{T}) where T <: BinaryElem = T(1)
Base.transpose(x::T) where T <: BinaryElem = x
Base.adjoint(x::T) where T <: BinaryElem = x
Base.conj(x::T) where T <: BinaryElem = x
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{T}) where T <: BinaryElem = T(rand(rng, poly_type(T)))

macro define_binary_elem(uint_size)
    gf2_elem_type = Symbol("BinaryElem$(uint_size)")
    poly_type = Symbol("BinaryPoly$(uint_size)")

    return quote
        struct $(gf2_elem_type) <: BinaryElem
            poly::$(poly_type)
        end
    end
end

@define_binary_elem 8
@define_binary_elem 16
@define_binary_elem 32
@define_binary_elem 64
@define_binary_elem 128

irreducible_poly(::Type{BinaryElem16}) = BinaryPoly16(UInt16(0b101101)) # a^16 + a^5 + a^3 + a^2 + 1, standard
irreducible_poly(::Type{BinaryElem32}) = BinaryPoly32(UInt32(0b11001 | 1 << 7 | 1 << 9 | 1 << 15)) # a^32 + a^15 + a^9 + a^7 + x^4 + x^3 + 1, Conway
irreducible_poly(::Type{BinaryElem128}) = BinaryPoly128(UInt128(0b10000111)) # x^128 + x^7 + x^2 + x + 1, standard (AES)

get_elem_type(::Type{BinaryPoly16}) = BinaryElem16
get_elem_type(::Type{BinaryPoly32}) = BinaryElem32
get_elem_type(::Type{BinaryPoly64}) = BinaryElem64
get_elem_type(::Type{BinaryPoly128}) = BinaryElem128

*(a::T, b::T) where T <: BinaryElem = mod_irreducible(poly(a)*poly(b))
+(a::T, b::T) where T <: BinaryElem = T(poly(a)+poly(b))

# when we compute irreducible / a then we have to be careful since irreducible requires 1 more bit than maximal bitsize of our polynomial
function div_irreducible(a::T) where {T<:BinaryElem}
    @assert binary_val(a) != 0

    shift = leading_zeros(binary_val(a)) + 1
    q0 = one(poly(a)) << shift
    r0 = irreducible_poly(T) + (poly(a) << shift)

    q, r = divrem(r0, poly(a))
    return q0 + q, r
end

function inv(a::T) where {T<:BinaryElem}
    q, r = div_irreducible(a) # p / a :: p = q*a + r 
    _, t, s = egcd(poly(a), r) # t*a + s*r = 1 = t*a + s*(p-q*a)
    # => (t-s*q)*a + s*p = 1
    # => ^^^^^^^ is inv of a modulo p, returned below

    T(t) + mod_irreducible(s * q)
end

@generated function compute_tmp(hi::T) where T <: BinaryPoly
    exprs = Expr[:(tmp = hi)]

    irr_poly = irreducible_poly(get_elem_type(T))
    
    for idx in 0:(bitsize(T)-1)
        if (binary_val(irr_poly) >> idx) & 1 == 1
            push!(exprs, :(tmp += hi >> $(bitsize(T) - idx)))
        end
    end

    push!(exprs, :(return tmp))

    return Expr(:block, exprs...)
end

@generated function compute_res(lo::T, tmp::T) where T <: BinaryPoly
    exprs = Expr[:(res = lo)]

    irr_poly = irreducible_poly(get_elem_type(T))
    
    for idx in 0:(bitsize(T)-1)
        if (binary_val(irr_poly) >> idx) & 1 == 1
            push!(exprs, :(res += tmp << $(idx)))
        end
    end

    push!(exprs, :(return res))

    return Expr(:block, exprs...)
end

function mod_irreducible(a::T) where T <: BinaryPoly
    (hi, lo) = split(a)

    tmp = compute_tmp(hi)
    res = compute_res(lo, tmp)

    Th = half_type(T)

    return get_elem_type(Th)(res)
end

function ^(a::T, n::U) where {T <: BinaryElem, U <: Integer}
    if a == zero(T)
        return zero(T)
    end

    t = one(T)

    while n > 0
        if n & 1 == 1
            t *= a
        end
        a *= a
        n >>= 1
    end

    return t
end


# beta is given as a root of F_2^m irreducible polynomial in F_2^128
# : beta = Fm.polynomial().roots(Fn, multiplicities=False)[0]
@generated function betas(::Type{T}) where T <: Union{BinaryElem16, BinaryElem32}
    bsize = bitsize(T)
    if T == BinaryElem16
        beta = BinaryElem128(UInt128(44320122245670141922313918651005395719))
    elseif T == BinaryElem32
        beta = BinaryElem128(UInt128(23246355947528323030879441634950214446))
    end
    bs = Vector{BinaryElem128}(undef, bsize)
    bs[1] = BinaryElem128(1)

    for (i, _) in enumerate(bs)
        if i == 1
            continue
        end
        bs[i] = bs[i-1]*beta
    end

    :($(bs))
end

promote_rule(::Type{BinaryElem128}, ::Type{T}) where T <: BinaryElem = BinaryElem128

function convert(::Type{BinaryElem128}, x::T) where T <: Union{BinaryElem16, BinaryElem32}
    result = BinaryElem128(0)
    bs = betas(T)
    x_bits = binary_val(x)

    for i in 1:bitsize(T)
        if (x_bits >> (i-1)) & 1 == 1
            result += bs[i]
        end
    end

    return result
end

# Note: The following code is unwrapped generic implementation of mod_irreducible
# function mod_irreducible(a::BinaryPoly64)
# irreducible = a^32 + a^15 + a^9 + a^7 + x^4 + x^3 + 1, standard
#     (hi, lo) = split(a)

#     tmp = hi + (hi >> (32 - 15)) + (hi >> (32 - 9)) + (hi >> (32 - 7)) + (hi >> (32 - 4)) + (hi >> (32 - 3))
#     res = lo + tmp + (tmp << 3) + (tmp << 4) + (tmp << 7) + (tmp << 9) + (tmp << 15)

#     return BinaryElem32(res)
# end