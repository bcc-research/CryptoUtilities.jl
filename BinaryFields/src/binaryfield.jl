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
    _, t, s = egcd(poly(a), r) # t*a + s*r = 1 = t*a + s*(p-q*a)
    # => (t-s*q)*a + s*p = 1
    # => ^^^^^^^ is inv of a modulo p, returned below

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

function xor_selected(x::UInt16, a::Vector{UInt128})
    result = UInt128(0)
    for i in 0:15
        if (x >> i) & 0x1 == 1
            result ⊻= a[i + 1]
        end
    end
    return result
end

function just_wrap(xs::Vector{UInt16}, as::Vector{Vector{UInt128}})
    n = length(xs)
    results = Vector{UInt128}(undef, n)
    for i in 1:n
        results[i] = xor_selected(xs[i], as[i])
    end
    return results
end
