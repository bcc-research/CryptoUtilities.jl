export BinaryPoly8, BinaryPoly16, BinaryPoly32, BinaryPoly64, BinaryPoly128

abstract type BinaryPoly end

# Generic stuff
binary_val(x::T) where T <: BinaryPoly = x.value
Base.zero(::T) where T <: BinaryPoly = T(0)
Base.zero(::Type{T}) where T <: BinaryPoly = T(0)
Base.transpose(x::T) where T <: BinaryPoly = x
Base.adjoint(x::T) where T <: BinaryPoly = x

macro define_binary_poly(uint_size)
    gf2_elem_type = Symbol("BinaryPoly$(uint_size)")
    uint_type = Symbol("UInt$(uint_size)")

    return quote
        struct $(gf2_elem_type) <: BinaryPoly
            value::$(uint_type)
        end

        $(esc(:primitive_type))(::Type{$(gf2_elem_type)}) = $(uint_type)
    end

end

@define_binary_poly 8
@define_binary_poly 16
@define_binary_poly 32
@define_binary_poly 64
@define_binary_poly 128

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{T}) where T <: BinaryPoly = T(rand(rng, primitive_type(T)))

Base.convert(::Type{T}, v::U) where {T <: BinaryPoly, U <: BinaryPoly} = T(binary_val(v))

+(a::T, b::T) where T <: BinaryPoly = T(binary_val(a) ⊻ binary_val(b))

function *(a::BinaryPoly64, b::BinaryPoly64)
    pmull_res = Vec(ccall("llvm.aarch64.neon.pmull64",
                          llvmcall,
                          NTuple{16, VecElement{UInt8}},
                          (UInt64, UInt64),
                          binary_val(a), binary_val(b)))

    return BinaryPoly128(reinterpret(UInt128, pmull_res))
end

# function split_long(a::UInt128)
#     UInt64(a >> 64), UInt64(a & typemax(UInt64))
# end

unwrap_value(a::VecElement{T}) where T = a.value

split(a::BinaryPoly128) = BinaryPoly64.(unwrap_value.(reinterpret(NTuple{2, VecElement{UInt64}}, binary_val(a))))
shift_upper_bits(a::BinaryPoly64) = BinaryPoly128(UInt128(binary_val(a)) << 64)

function *(a::BinaryPoly128, b::BinaryPoly128)
    a_lo, a_hi = split(a)
    b_lo, b_hi = split(b)

    z0 = a_lo * b_lo
    z1 = (a_lo + a_hi) * (b_lo + b_hi)
    z2 = a_hi * b_hi

    result_lo = z0
    result_hi = z2
    result_mid_lo, result_mid_hi = split(z0 + z1 + z2)

    lo_bits = shift_upper_bits(result_mid_lo) + result_lo
    hi_bits = result_hi + convert(BinaryPoly128, result_mid_hi)

    # XXX: Oh man, this is different from the way `reinterpret` works in Julia.
    return (hi_bits, lo_bits)
end

# Can make generic, fine for now, see docs for why this works
function reduce_step(a::NTuple{2, BinaryPoly128}, poly)
    a_hi, a_lo = a
    res_hi, res_lo = a_hi * poly

    return (res_hi, a_lo + res_lo)
end

function reduce_poly(a::NTuple{2, BinaryPoly128}, poly)
    a = reduce_step(a, poly)
    _, a_lo = reduce_step(a, poly)
    return a_lo
end

# # XXX: Maybe should be a different type, but good enough for now.  Computes the
# # divisor and remainder of a / b, interpreting the bits as coefficients of a
# # polyomial over F_2. Used for inversion. Also assumes that (for now) things are
# # not silly i.e., that `a` isn't mostly zeros. Can be optimized in that case.
# function divrempoly(a::UInt128, b::UInt128)
#     @assert b != 0

#     shift = leading_zeros(b)
#     q = UInt128(0)

#     bit_post = UInt128(1) << UInt128(127)
#     bit_post_div = UInt128(1) << shift
#     b <<= shift

#     while shift >= 0
#         if (a & bit_post) != 0
#             q |= bit_post_div
#             a ⊻= b
#         end
#         shift -= 1
#         b >>= 1
#         bit_post >>= 1
#         bit_post_div >>= 1
#     end

#     return q, a
# end

# # when we compute irreducible / a then we have to be careful since irreducible requires more than 128 bits
# function div_irreducible(a::UInt128, irr_low)
#     @assert a != 0 

#     shift = leading_zeros(a) + 1 
#     q0 = UInt128(1) << shift
#     r0 = irr_low ⊻ (a << shift)

#     (q, r) = divrempoly(r0, a)

#     return (q0 ⊻ q, r)
# end

# # g, t, s = egcd(a, b) => t * a + s * b = g = gcd(a, b)
# # | via recursion: t' * b + s' * (a % b) = gcd(a, b)
# # | => t' * b + s' * (a + a/b * b)
# # | => s' * a + t' + (s' * a/b) * b = gcd(a, b)
# function egcd(r_1::UInt128, r_2::UInt128)
#     if r_2 == 0
#         @assert r_1 != 0
#         return r_1, UInt128(1), UInt128(0)
#     else
#         q, r_3 = divrempoly(r_1, r_2)
#         g, t, s = egcd(r_2, r_3)
#         _, qs = carryless_mul(q, s)
#         return g, s, qs ⊻ t
#         g, x1, y1 = egcd(r_2, r_3)
#         _, mul_res = carryless_mul(q, y1)
#         return g, y1, mul_res ⊻ x1
#     end
# end

# # XXX: We should make these generic
# function inv(a::GF2_128Elem)
#     q, r = div_irreducible(binary_val(a), irreducible_poly(a)) # p / a :: p = q*a + r 
#     _, t, s = egcd(binary_val(a), r)
#     # t*a + s*r = 1 = t*a + s*(p-q*a)
#     # => (t-s*q)*a + s*p = 1
#     # => ^^^^^^^ is inv of a
#     return GF2_128Elem(t) + GF2_128Elem(s)*GF2_128Elem(q)
# end

# # --- Other elements of other sizes ---
# macro define_GF2_Elem(uint_size)
#     gf2_elem_type = Symbol("GF2_$(uint_size)Elem")
#     uint_type = Symbol("UInt$(uint_size)")

#     return quote
#         struct $(gf2_elem_type) <: BinaryFieldElem
#             value::$(uint_type)
#         end
#     end
# end

# export GF2_8Elem
# # Right now, just handle these by upconverting to GF2_128Elem
# @define_GF2_Elem 8
# @define_GF2_Elem 16
# @define_GF2_Elem 32
# @define_GF2_Elem 64

# # XXX: Make this fast via repeated squaring
# function Base.convert(::Type{GF2_128Elem}, v::T) where T <: BinaryFieldElem
#     a = v.value
#     # Should be constant, check decompilation
#     bit_size = 8*sizeof(v)
#     skip_size = div(128, bit_size)

#     output = UInt128(0)
#     for i in 1:bit_size
#         curr_idx = (a >> (i-1)) & 1
#         output |= UInt128(curr_idx) << UInt128(skip_size*(i-1))
#     end

#     GF2_128Elem(output)
# end
