
function one_step_fft(v::Vector{GF2_128Elem}; twiddles) # v is of length 4, GF2_128
    i0 = FFTMatrix(twiddles[1])*v
    
    i10 = FFTMatrix(twiddles[2])*i0[1:2]
    i11 = FFTMatrix(twiddles[3])*i0[3:4]

    i0
end

export FFTMatrix

# 2length(D) x 2length(D)
# [I   λI  ]
# [I (λ+1)I]
struct FFTMatrix3{T}
    λ::T
end

# TODO: Method to pretty print the matrix

FFTMatrix = FFTMatrix3

function *(M::FFTMatrix{T}, v::Vector{T}) where T
    n = length(v)
    @assert n % 2 == 0
    n_div2 = div(n, 2)
    @views u, w = v[1:n_div2], v[n_div2+1:end]

    udw = u + M.λ*w

    return [udw; udw + w]
end