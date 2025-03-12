# Using SIMD for fast reinterpret. Julia's base `reinterpret` is insanely slow
# in this particular case likely due to (in this case unnecessary) type safety
# checks

abstract type BinaryElem <: Number end

# --- Other elements of other sizes ---
# macro define_binary_elem(uint_size)
#     gf2_elem_type = Symbol("BinaryElem$(uint_size)")
#     uint_type = Symbol("UInt$(uint_size)")

#     return quote
#         struct $(gf2_elem_type) <: BinaryElem
#             value::$(uint_type)
#         end
#     end
# end

# @define_binary_elem 8
# @define_binary_elem 16
# @define_binary_elem 32
# @define_binary_elem 64

poly(x::T) where T <: BinaryElem = x.poly
binary_val(x::T) where T <: BinaryElem = binary_val(poly(x))
Base.zero(::T) where T <: BinaryElem = T(0)
Base.zero(::Type{T}) where T <: BinaryElem = T(0)
Base.transpose(x::T) where T <: BinaryElem = x
Base.adjoint(x::T) where T <: BinaryElem = x

Base.convert(::Type{T}, v::T) where {T<:BinaryElem} = v

*(a::T, b::T) where T<: BinaryElem = mod_irreducible(poly(a)*poly(b))
+(a::T, b::T) where T<: BinaryElem = T(poly(a)+poly(b))
<<(a::T, n::Int) where {T<:BinaryElem} = T(poly(a) << n)
>>(a::T, n::Int) where {T<:BinaryElem} = T(poly(a) >> n)

# when we compute irreducible / a then we have to be careful since irreducible requires 1 more bit than maximal bitsize of our polynomial
function div_irreducible(a::T) where {T<:BinaryElem}
    @assert binary_val(a) != 0

    shift = leading_zeros(binary_val(a)) + 1
    q0 = one(T) << shift
    r0 = irreducible_poly(T) + poly((a << shift))

    q, r = divrem(r0, poly(a))
    return q0 + q, r
end

function inv(a::T) where {T<:BinaryElem}
    q, r = div_irreducible(a) # p / a :: p = q*a + r 
    _, t, s = egcd(poly(a), r)
    # t*a + s*r = 1 = t*a + s*(p-q*a)
    # => (t-s*q)*a + s*p = 1
    # => ^^^^^^^ is inv of a
    convert(T, t) + mod_irreducible(s * q)
end
