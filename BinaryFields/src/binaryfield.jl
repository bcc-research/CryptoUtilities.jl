# Using SIMD for fast reinterpret. Julia's base `reinterpret` is insanely slow
# in this particular case likely due to (in this case unnecessary) type safety
# checks

abstract type BinaryElem <: Number end

poly(x::T) where T <: BinaryElem = x.poly
binary_val(x::T) where T <: BinaryElem = binary_val(poly(x))
Base.zero(::T) where T <: BinaryElem = T(0)
Base.zero(::Type{T}) where T <: BinaryElem = T(0)
Base.transpose(x::T) where T <: BinaryElem = x
Base.adjoint(x::T) where T <: BinaryElem = x

Base.convert(::Type{T}, v::T) where {T<:BinaryElem} = v

# Define a GF2^128 element of a binary field
struct BinaryElem128 <: BinaryElem
    poly::BinaryPoly128
end
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{BinaryElem128}) = BinaryElem128(rand(rng, BinaryPoly128))

irreducible_poly(::BinaryElem128) = UInt128(0b10000111) # x^128 + x^7 + x^2 + x + 1, standard

# Make these macro-generated
function mod_irreducible(a::NTuple{2, BinaryPoly128})
    (hi, lo) = a

    tmp = hi + (hi >> 127) + (hi >> 126) + (hi >> 121)
    res = lo + tmp + (tmp << 1) + (tmp << 2) + (tmp << 7)

    return res
end

*(a::BinaryElem128, b::BinaryElem128) = mod_irreducible(poly(a)*poly(b))
+(a::T, b::T) where T <: BinaryElem = BinaryElem128(poly(a)+poly(b))

# XXX: Maybe should be a different type, but good enough for now.  Computes the
# divisor and remainder of a / b, interpreting the bits as coefficients of a
# polyomial over F_2. Used for inversion. Also assumes that (for now) things are
# not silly i.e., that `a` isn't mostly zeros. Can be optimized in that case.
function divrempoly(a::UInt128, b::UInt128)
    @assert b != 0

    shift = leading_zeros(b)
    q = UInt128(0)

    bit_post = UInt128(1) << UInt128(127)
    bit_post_div = UInt128(1) << shift
    b <<= shift

    while shift >= 0
        if (a & bit_post) != 0
            q |= bit_post_div
            a ⊻= b
        end
        shift -= 1
        b >>= 1
        bit_post >>= 1
        bit_post_div >>= 1
    end

    return q, a
end

# when we compute irreducible / a then we have to be careful since irreducible requires more than 128 bits
function div_irreducible(a::UInt128, irr_low)
    @assert a != 0 

    shift = leading_zeros(a) + 1 
    q0 = UInt128(1) << shift
    r0 = irr_low ⊻ (a << shift)

    (q, r) = divrempoly(r0, a)

    return (q0 ⊻ q, r)
end

# g, t, s = egcd(a, b) => t * a + s * b = g = gcd(a, b)
# | via recursion: t' * b + s' * (a % b) = gcd(a, b)
# | => t' * b + s' * (a + a/b * b)
# | => s' * a + t' + (s' * a/b) * b = gcd(a, b)
function egcd(r_1::UInt128, r_2::UInt128)
    if r_2 == 0
        @assert r_1 != 0
        return r_1, UInt128(1), UInt128(0)
    else
        q, r_3 = divrempoly(r_1, r_2)
        g, t, s = egcd(r_2, r_3)
        _, qs = carryless_mul(q, s)
        return g, s, qs ⊻ t
        g, x1, y1 = egcd(r_2, r_3)
        _, mul_res = carryless_mul(q, y1)
        return g, y1, mul_res ⊻ x1
    end
end

# XXX: We should make these generic
function inv(a::BinaryElem128)
    q, r = div_irreducible(binary_val(a), irreducible_poly(a)) # p / a :: p = q*a + r 
    _, t, s = egcd(binary_val(a), r)
    # t*a + s*r = 1 = t*a + s*(p-q*a)
    # => (t-s*q)*a + s*p = 1
    # => ^^^^^^^ is inv of a
    return BinaryElem128(t) + BinaryElem128(s)*BinaryElem128(q)
end

# --- Other elements of other sizes ---
macro define_binary_elem(uint_size)
    gf2_elem_type = Symbol("BinaryElem$(uint_size)")
    uint_type = Symbol("UInt$(uint_size)")

    return quote
        struct $(gf2_elem_type) <: BinaryElem
            value::$(uint_type)
        end
    end
end

# Right now, just handle these by upconverting to BinaryElem128
@define_binary_elem 8
@define_binary_elem 16
@define_binary_elem 32
@define_binary_elem 64

# XXX: Make this fast via repeated squaring
# function Base.convert(::Type{GF2_128Elem}, v::T) where T <: BinaryFieldElem
#     a = v.value
#     # Should be constant, check decompilation
#     @show bit_size = 8*sizeof(v)
#     skip_size = div(128, bit_size)

#     output = UInt128(0)
#     for i in 1:bit_size
#         curr_idx = (a >> (i-1)) & 1
#         output |= UInt128(curr_idx) << UInt128(skip_size*(i-1))
#     end

#     GF2_128Elem(output)
# end