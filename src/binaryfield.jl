import Random

# Using SIMD for fast reinterpret. Julia's base `reinterpret` is insanely slow
# in this particular case likely due to (in this case unnecessary) type safety
# checks
using SIMD
import Base: *, +

abstract type BinaryFieldElem <: Number end

binary_val(x::T) where T <: BinaryFieldElem = x.value
Base.zero(::T) where T <: BinaryFieldElem = T(0)
Base.zero(::Type{T}) where T <: BinaryFieldElem = T(0)
Base.transpose(x::T) where T <: BinaryFieldElem = x
Base.adjoint(x::T) where T <: BinaryFieldElem = x

export GF2_128Elem

# Define a GF2^128 element of a binary field
struct GF2_128Elem <: BinaryFieldElem
    value::UInt128
end
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{GF2_128Elem}) = GF2_128Elem(rand(rng, UInt128))


irreducible_poly(::GF2_128Elem) = UInt128(0b10000111) # x^128 + x^7 + x^2 + x + 1, standard

# Only works for ARM chips with crypto NEON extensions
function carryless_mul(a::UInt64, b::UInt64)
    pmull_res = Vec(ccall("llvm.aarch64.neon.pmull64",
                          llvmcall,
                          NTuple{16, VecElement{UInt8}},
                          (UInt64,UInt64),
                          a, b))

    return reinterpret(UInt128, pmull_res)
end

function split_long(a::UInt128)
    UInt64(a >> 64), UInt64(a & typemax(UInt64))
end

function carryless_mul(a::UInt128, b::UInt128)
    a_hi, a_lo = split_long(a)
    b_hi, b_lo = split_long(b)

    z0 = carryless_mul(a_lo, b_lo)
    z1 = carryless_mul(a_lo ⊻ a_hi, b_lo ⊻ b_hi)
    z2 = carryless_mul(a_hi, b_hi)

    result_lo = z0
    result_hi = z2
    result_mid = z0 ⊻ z1 ⊻ z2

    lo_bits = (result_mid << 64) ⊻ result_lo
    hi_bits = result_hi ⊻ (result_mid >> 64)
    
    return (hi_bits, lo_bits)
end

# Can make generic, fine for now
function reduce_once(a::NTuple{2, UInt128}, poly)
    a_hi, a_lo = a
    res_hi, res_lo = carryless_mul(a_hi, poly)

    return (res_hi, a_lo ⊻ res_lo)
end

function reduce_poly(a::NTuple{2, UInt128}, poly)
    a = reduce_once(a, poly)
    _, a_lo = reduce_once(a, poly)
    return a_lo
end

function *(a::T, b::GF2_128Elem) where T <: BinaryFieldElem
    a_upconv = convert(GF2_128Elem, a)
    c_mul = carryless_mul(binary_val(a_upconv), binary_val(b))
    return GF2_128Elem(reduce_poly(c_mul, irreducible_poly(b)))
end

function +(a::GF2_128Elem, b::GF2_128Elem)
    GF2_128Elem(a.value ⊻ b.value)
end

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
function div_irreducible(a::UInt128)
    @assert a != 0 
    irr_low = UInt128(0b10000111) # x^7 + x^2 + x + 1

    shift = leading_zeros(a) + 1 
    q0 = UInt128(1) << shift
    r0 = irr_low ⊻ (a << shift)

    (q, r) = divrempoly(r0, a)

    return (q0 ⊻ q, r)
end

# egcd(a, b) => t * a + s * b = gcd(a, b)
# | via recursion: t' * b + s' * (a % b) = gcd(a, b)
# | => t' * b + s' * (a + a/b * b)
# | => s' * a + t' + (s' * a/b) * b = gcd(a, b)
function egcd(r_1::UInt128, r_2::UInt128)
    if r_2 == UInt128(0)
        @assert r_1 != 0
        return r_1, UInt128(1), UInt128(0)
    else
        q, r_3 = divrempoly(r_1, r_2)
        g, t, s = egcd(r_2, r_3)
        _, qs = carryless_mul(q, s)
        return g, s, qs ⊻ t
    end
end


function inv(a::GF2_128Elem)
    q, r = div_irreducible(binary_val(a)) # p / a :: p = q*a + r 
    _, t, s = egcd(binary_val(a), r)
    return GF2_128Elem(t) + GF2_128Elem(s)*GF2_128Elem(q)
end

# --- Other elements of other sizes ---
macro define_GF2_Elem(uint_size)
    gf2_elem_type = Symbol("GF2_$(uint_size)Elem")
    uint_type = Symbol("UInt$(uint_size)")

    return quote
        struct $(gf2_elem_type) <: BinaryFieldElem
            value :: $(uint_type)
        end
    end
end

export GF2_8Elem
# Right now, just handle these by upconverting to GF2_128Elem
@define_GF2_Elem 8
@define_GF2_Elem 16
@define_GF2_Elem 32
@define_GF2_Elem 64

# XXX: Make this fast via repeated squaring
function Base.convert(::Type{GF2_128Elem}, v::T) where T <: BinaryFieldElem
    a = v.value
    # Should be constant, check decompilation
    bit_size = 8*sizeof(v)
    skip_size = div(128, bit_size)

    output = UInt128(0)
    for i in 1:bit_size
        curr_idx = (a >> (i-1)) & 1
        output |= UInt128(curr_idx) << UInt128(skip_size*(i-1))
    end

    GF2_128Elem(output)
end
