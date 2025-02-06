using Base.Threads

is_pow_2(n) = 2^(round(Int, log2(n))) == n

function fft_twiddles!(v; twiddles, idx=1)
    if length(v) == 1
        return v
    end

    mul_inplace!(v, twiddles[idx])
    
    u, v = split_half(v)

    @views fft_twiddles!(u; twiddles, idx=2*idx)
    @views fft_twiddles!(v; twiddles, idx=2*idx+1)
end

function fft_twiddles_parallel!(v; twiddles, idx=1, thread_depth=nothing)
    if length(v) == 1
        return v
    end

    if isnothing(thread_depth)
        thread_depth = round(Int, log2(nthreads()))
        @info "Setting thread depth to $thread_depth"
    end

    mul_inplace!(v, twiddles[idx])
    
    u, v = split_half(v)

    if thread_depth > 0
        @threads for j in 1:2
            @views fft_twiddles_parallel!(j==1 ? u : v; twiddles, idx=2*idx+j-1, thread_depth=thread_depth-1)
        end
    else 
        @views fft_twiddles!(u; twiddles, idx=2*idx)
        @views fft_twiddles!(v; twiddles, idx=2*idx+1)
    end
end

# Takes a vector `v` of even length and returns the
# first and second halves (as views)
function split_half(v)
    n = length(v)
    n_div2 = div(n, 2)
    return @views v[1:n_div2], v[n_div2+1:end]
end

# Multiplies (in place!) the FFT matrix with twiddle
# factor λ and vector `v`
# [I   λI  ]
# [I (λ+1)I]
function mul_inplace!(v, λ)
    u, w = split_half(v)
    @views begin
        @. u += λ*w
        w .+= u
    end
end