is_pow_2(n) = 2^(round(Int, log2(n))) == n

function fft_twiddles!(v; twiddles, idx=1)
    if length(v) == 1
        return v
    end

    mul_inplace!(v, twiddles[idx])
    
    u, w = split_half(v)

    @views fft_twiddles!(u; twiddles=twiddles, idx=2*idx)
    @views fft_twiddles!(w; twiddles=twiddles, idx=2*idx+1)
end

function split_half(v)
    n = length(v)
    n_div2 = div(n, 2)
    return @views v[1:n_div2], v[n_div2+1:end]
end

# [I   λI  ]
# [I (λ+1)I]
function mul_inplace!(v, λ)
    u, w = split_half(v)
    @views begin
        @. u += λ*w
        w .+= u
    end
end