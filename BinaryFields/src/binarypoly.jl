
abstract type BinaryPoly end

# Generic stuff
Base.zero(::T) where {T<:BinaryPoly} = T(0)
Base.zero(::Type{T}) where {T<:BinaryPoly} = T(0)
Base.one(::T) where {T<:BinaryPoly} = T(1)
Base.one(::Type{T}) where {T<:BinaryPoly} = T(1)
Base.transpose(x::T) where {T<:BinaryPoly} = x
Base.adjoint(x::T) where {T<:BinaryPoly} = x

macro define_binary_poly(uint_size)
    gf2_elem_type = Symbol("BinaryPoly$(uint_size)")
    uint_type = Symbol("UInt$(uint_size)")

    return quote
        struct $(gf2_elem_type) <: BinaryPoly
            value::$(uint_type)
        end
    end

end

@define_binary_poly 8
@define_binary_poly 16
@define_binary_poly 32
@define_binary_poly 64
@define_binary_poly 128

struct BinaryPoly256 <: BinaryPoly
    value::Tuple{UInt128, UInt128}
end

# Stuff that depends on struct definitions
binary_val(x::T) where {T<:BinaryPoly} = x.value
primitive_type(::Type{T}) where T <: BinaryPoly = fieldtypes(T)[1]

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{T}) where {T<:BinaryPoly} = T(rand(rng, primitive_type(T)))

Base.convert(::Type{T}, v::U) where {T<:BinaryPoly,U<:BinaryPoly} = T(binary_val(v))
Base.convert(::Type{BinaryPoly256}, v::NTuple{2, BinaryPoly128}) = T(binary_val.(v))
Base.convert(::Type{T}, x) where {T<:BinaryPoly} = T(x)

+(a::T, b::T) where {T<:BinaryPoly} = T(binary_val(a) ⊻ binary_val(b))

<<(a::T, n::Int) where {T<:BinaryPoly} = T(binary_val(a) << n)
>>(a::T, n::Int) where {T<:BinaryPoly} = T(binary_val(a) >> n)

saturate(::Type{T}, a::U) where {T, U <: BinaryPoly} = T(binary_val(a) & typemax(primitive_type(T)))
function join(a::T, b::T) where {T <: BinaryPoly} 
    T_double = double_type(T)

    a_double = convert(T_double, a)
    a_double <<= sizeof(primitive_type(T)) * 8
    b_double = convert(T_double, b)

    return a_double + b_double
end

double_type(::Type{BinaryPoly8}) = BinaryPoly16
double_type(::Type{BinaryPoly16}) = BinaryPoly32
double_type(::Type{BinaryPoly32}) = BinaryPoly64
double_type(::Type{BinaryPoly64}) = BinaryPoly128

half_type(::Type{BinaryPoly16}) = BinaryPoly8
half_type(::Type{BinaryPoly32}) = BinaryPoly16
half_type(::Type{BinaryPoly64}) = BinaryPoly32
half_type(::Type{BinaryPoly128}) = BinaryPoly64

function split(x::T) where {T <: BinaryPoly} 
    T_half = half_type(T)

    lo = saturate(T_half, x)
    hi = convert(T_half, x >> (sizeof(T_half) * 8))

    hi, lo
end

split(x::BinaryPoly256) = BinaryPoly128.(binary_val(x))

shift_upper_bits(a::BinaryPoly64) = BinaryPoly128(UInt128(binary_val(a)) << 64)

@generated function *(a::BinaryPoly64, b::BinaryPoly64)
    if Sys.ARCH == :aarch64
        quote
            pmull_res = Vec(ccall("llvm.aarch64.neon.pmull64",
                llvmcall,
                NTuple{16,VecElement{UInt8}},
                (UInt64, UInt64),
                binary_val(a), binary_val(b)
            ))

            return BinaryPoly128(reinterpret(UInt128, pmull_res))
        end
    elseif Sys.ARCH == :x86_64
        # XXX: Might need some optimizations, not sure these are the best types
        # right now!
        # XXX: Also, should rewrite * depending on arch? Note that we can
        # multiply any parts of the register using the last argument to
        # pclmulqdq, which should reduce data movement.
        quote
            pmull_res = Vec(ccall("llvm.x86.pclmulqdq",
                llvmcall,
                NTuple{2,VecElement{UInt64}},
                (NTuple{2,VecElement{UInt64}}, NTuple{2,VecElement{UInt64}}, UInt8),
                VecElement{UInt64}.((binary_val(a), 0)), VecElement{UInt64}.((binary_val(b), 0)), 0x00
            ))

            return BinaryPoly128(reinterpret(UInt128, pmull_res))
        end
    else
        error("Unsupported architecture for binary polynomial 64-bit multiplication")
    end
end

function *(a::BinaryPoly128, b::BinaryPoly128)
    a_hi, a_lo = split(a)
    b_hi, b_lo = split(b)

    z0 = a_lo * b_lo
    z1 = (a_lo + a_hi) * (b_lo + b_hi)
    z2 = a_hi * b_hi

    result_lo = z0
    result_hi = z2
    result_mid_hi, result_mid_lo = split(z0 + z1 + z2)

    lo_bits = shift_upper_bits(result_mid_lo) + result_lo
    hi_bits = result_hi + convert(BinaryPoly128, result_mid_hi)

    return BinaryPoly256(binary_val.((hi_bits, lo_bits)))
end

Poly2x16 = NTuple{2, BinaryPoly16}

function embed_p64(a::Poly2x16)
    # This should compile down to two loads on a quad register on aarch64?
    res = UInt64(binary_val(a[2])) << 32 | binary_val(a[1])

    return BinaryPoly64(res)
end

"""
Perform 2x16 multiplication of a, b via 64 bit multiplication using
(0 | a[2] | 0 | a[1] ) * (0 | b[2] | 0 | b[1]) and then unpacking
"""
function *(a::Poly2x16, b::Poly2x16)
    a_64, b_64 = embed_p64(a), embed_p64(b)
    res = a_64*b_64

    res1 = saturate(BinaryPoly32, res)
    res2 = saturate(BinaryPoly32, res >> 64)

    return (res1, res2)
end

function *(λ::BinaryPoly16, b::Poly2x16)
    λ_a = (λ, λ)

    return λ_a * b
end

function *(a::T, b::T) where {T <: Union{BinaryPoly8, BinaryPoly16, BinaryPoly32}}
    res = convert(BinaryPoly64, a) * convert(BinaryPoly64, b)

    return convert(double_type(T), res)
end

# Computes the divisor and remainder of a / b
function divrem(a::T, b::T) where {T <: BinaryPoly}
    @assert binary_val(b) != 0

    shift = leading_zeros(binary_val(b))
    q = T(0)

    bit_post = T(1) << (sizeof(primitive_type(T)) * 8 - 1)
    bit_post_div = T(1) << shift
    b = b << shift

    while shift >= 0
        if (binary_val(a) & binary_val(bit_post)) != 0
            q = q + bit_post_div
            a = a + b
        end
        shift -= 1

        b >>= 1
        bit_post >>= 1
        bit_post_div >>= 1
    end

    return q, a
end

# g, t, s = egcd(a, b) => t * a + s * b = g = gcd(a, b)
# | via recursion: t' * b + s' * (a % b) = gcd(a, b)
# | => t' * b + s' * (a + a/b * b)
# | => s' * a + t' + (s' * a/b) * b = gcd(a, b)
function egcd(r_1::T, r_2::T) where {T<:BinaryPoly}
    if binary_val(r_2) == 0
        @assert binary_val(r_1) != 0
        return r_1, T(1), T(0)
    else
        q, r_3 = divrem(r_1, r_2)
        g, t, s = egcd(r_2, r_3)
        _, qs = split(q * s)
        return g, s, qs + t
    end
end
