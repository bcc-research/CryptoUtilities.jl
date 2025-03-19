# Using SIMD for fast reinterpret. Julia's base `reinterpret` is insanely slow
# in this particular case likely due to (in this case unnecessary) type safety
# checks

abstract type BinaryElem <: Number end

poly(x::T) where T <: BinaryElem = x.poly
poly_type(::Type{T}) where T <: BinaryElem = fieldtypes(T)[1]
binary_val(x::T) where T <: BinaryElem = binary_val(poly(x))
Base.zero(::T) where T <: BinaryElem = T(0)
Base.zero(::Type{T}) where T <: BinaryElem = T(0)
Base.one(::T) where T <: BinaryElem = T(1)
Base.one(::Type{T}) where T <: BinaryElem = T(1)
Base.transpose(x::T) where T <: BinaryElem = x
Base.adjoint(x::T) where T <: BinaryElem = x
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{T}) where T <: BinaryElem = T(rand(rng, poly_type(T)))

#tmp 
export binary_val

*(a::T, b::T) where T<: BinaryElem = mod_irreducible(poly(a)*poly(b))
+(a::T, b::T) where T<: BinaryElem = T(poly(a)+poly(b))

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
    _, t, s = egcd(poly(a), r)
    # t*a + s*r = 1 = t*a + s*(p-q*a)
    # => (t-s*q)*a + s*p = 1
    # => ^^^^^^^ is inv of a
    T(t) + mod_irreducible(s * q)
end

function embed(a::T, ::Type{U}) where {T<:BinaryElem, U<:BinaryElem}
    @assert sizeof(T) < sizeof(U)
    k = sizeof(U) ÷ sizeof(T)

    res = poly(a)
    for _ in 1:k 
        res *= res
    end

    U(res)
end