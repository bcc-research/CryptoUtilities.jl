import Random

# Using SIMD for fast reinterpret. Julia's base `reinterpret` is insanely slow
# in this particular case likely due to (in this case unnecessary) type safety
# checks
using SIMD
import Base: *, +

abstract type BinaryFieldElem <: Number end

# | I 0 | | I D |
# | I I | | 0 I |  *v
#
#
# | I   D   |
# | I D + I |

binary_val(x::T) where T <: BinaryFieldElem = x.value
Base.zero(::T) where T <: BinaryFieldElem = T(0)
Base.zero(::Type{T}) where T <: BinaryFieldElem = T(0)
Base.transpose(x::T) where T <: BinaryFieldElem = x
Base.adjoint(x::T) where T <: BinaryFieldElem = x

function *(a::Bool, b::T) where T <: BinaryFieldElem
    a ? b : zero(b)
end

function *(a::T, M::UniformScaling{Bool}) where T <: BinaryFieldElem
    M.λ ? UniformScaling(a) : UniformScaling(false)
end

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

# XXX: Maybe should be a different type, but good enough for now.
# Computes the divisor and remainder of a / b, interpreting
# the bits as coefficients of a polyomial over F_2. Used for
# inversion.
function divrempoly_(a::UInt128, b::UInt128)
    shift = leading_zeros(b)
    curr_idx = 128
    b <<= shift
    q, r = UInt128(0), UInt128(0)
    while shift >= 0
        if curr_idx 
        shift -= 1
        end
    end

    return q, r
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
