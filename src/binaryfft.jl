using Base.Threads

export fft!, ifft!

"""
    fft_twiddles!(v; twiddles, idx=1)

In-place recursive FFT step using pre-calculated twiddle factors.

Recursively applies twiddle factors to vector `v` in-place, modifying it directly.
`v` must be a power of 2 in length.
# Arguments
- `v`: Vector to transform in-place (must be power of 2).
- `twiddles`: Pre-calculated twiddle factors.
- `idx`: Twiddle factor index (defaults to 1).
# Keyword Arguments
- `twiddles`:  Required, the collection of twiddle factors.
- `idx`: Optional, the index into the twiddle factor collection, defaults to 1.

# Notes
- Building block for FFT implementations.
- Operates in-place on `v`.
- `idx` is used internally to navigate `twiddles` during recursion.
"""
function fft_twiddles!(v; twiddles, idx=1)
    if length(v) == 1
        return v
    end

    fft_mul!(v, twiddles[idx])
    
    u, v = split_half(v)

    @views fft_twiddles!(u; twiddles, idx=2*idx)
    @views fft_twiddles!(v; twiddles, idx=2*idx+1)
end

# TODO: Add nice comments here
function ifft_twiddles!(v; twiddles, idx=1)
    if length(v) == 1
        return v
    end

    h1, h2 = split_half(v)

    @views ifft_twiddles!(h1; twiddles, idx=2*idx)
    @views ifft_twiddles!(h2; twiddles, idx=2*idx+1)

    ifft_mul!(v, twiddles[idx])
end

"""
    fft_twiddles_parallel!(v; twiddles, idx=1, thread_depth=nothing)

Recursively apply twiddle factors in parallel to perform a Fast Fourier Transform (FFT) step in-place.

Parallel recursive FFT step using pre-calculated twiddle factors. Operates in-place on vector `v` using Julia threads.
# Arguments
- `v`: Vector to transform in-place (length must be power of 2).
- `twiddles`: Pre-calculated twiddle factors.
- `idx`: Twiddle factor index (defaults to 1).
- `thread_depth`: Depth of parallel recursion. Defaults to `log2(nthreads())` initially, decrements recursively, serial below 0.
# Notes
- Parallel version of `fft_twiddles!` for multi-core speedup.
- `thread_depth` controls parallelism; optimize to balance overhead and benefit.
- Default `thread_depth` is inferred and logged if not provided.
"""
function fft_twiddles_parallel!(v; twiddles, idx=1, thread_depth=nothing)
    if length(v) == 1
        return v
    end

    if isnothing(thread_depth)
        thread_depth = round(Int, log2(nthreads()))
        if thread_depth > 0
          @info "Setting thread depth to $thread_depth"
        else
          @info "Setting thread depth to 0 (did you launch julia with `--threads=[thread_count]`?)"
        end
    end

    fft_mul!(v, twiddles[idx])
    
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

function split_half(v)
    n = length(v)
    n_div2 = div(n, 2)
    return @views v[1:n_div2], v[n_div2+1:end]
end

"""
    fft_mul!(v, λ)

Multiply in-place the FFT matrix with twiddle factor `λ` and vector `v`.

This function performs the following matrix-vector multiplication in-place:

[u'] = [I   λI  ][u]
[w'] = [I (λ+1)I][w]

where `I` is an identity matrix of size `length(v) ÷ 2` and `v` is a vector
partitioned into two halves `u` and `w` such that `v = [u; w]`.  The
multiplication is performed in-place, modifying the input vector `v` directly.

This operation is a step in the binary field FFT algorithm of XXX.

# Arguments
- `v`: A vector to be multiplied in-place. Modified upon function return.
- `λ`: The twiddle factor, a scalar value used in the matrix multiplication.
"""
# TODO: Maybe we should rename this function to something more convenient
function fft_mul!(v, λ)
    u, w = split_half(v)
    @views begin
        @. u += λ*w
        w .+= u
    end
end

# Ok this kind of works now, but it is slower than the above for some reason.
# Probably for later.
function fft_mul!(v::T, λ::BinaryElem16) where {T <: AbstractVector{BinaryElem16}}
    u, w = split_half(v)
    if length(u) > 1
        @views begin
            BinaryFields.scale_inplace_simd!(u, λ, v)
            w .+= u
        end
    else
        @views begin
            @. u += λ*w
            w .+= u
        end
    end
end

# TODO: Add nice comments here
function ifft_mul!(v, λ)
    lo, hi = split_half(v)
    @views begin
        hi .+= lo
        @. lo += λ*hi
    end
end


next_s(s_prev, s_prev_at_root) = s_prev*s_prev + s_prev_at_root*s_prev

# s_i(v_i) = (s_{i-1}(v_i))^2 - s_{i-1}(v_{i - 1})*s_{i-1}(v_i)
# note that first two elements of each layer are: 
# s_i(beta), s_i(beta + v_{i + 1}) 
# thus we can compute s_{i+1}(v_{i + 1}) as: 
# Example: say we want to compute s2(v2)
# s2(v2) = s1(v2)^2 - s1(v1)*s1(v2)
# prev layer starts with: [s1(b), s1(b + v2), ...] 
# so we can compute: s1(v2) = s1(b + v2) - s1(b) = s1(b) + s1(v2) - s1(b)
compute_s_at_root(prev_layer, s_prev_at_root) = next_s(prev_layer[2] + prev_layer[1], s_prev_at_root)

function layer_0!(layer, beta, k)
    T = eltype(layer)
    for i in 1:2^(k-1)
      l0i = beta
      l0i += T((i-1) << 1)
      layer[i] = l0i
    end

    # s0(v0)
    return T(1)
end

function layer_i!(layer, layer_len, s_prev_at_root)
    prev_layer_len = 2 * layer_len
    s_at_root = compute_s_at_root(layer, s_prev_at_root)
    for (idx, s_prev) in enumerate(@views layer[1:2:prev_layer_len])
      layer[idx] = next_s(s_prev, s_prev_at_root)
    end
    return s_at_root
end 

# TODO: add comments how pi polynomials have a similar subspace structure as full F_2^m space 
function precompute_basis_scalers_at_layer_k!(p, k, sk_at_vk)
	current_len = 2^k 
	p[current_len+1:2*current_len] .= sk_at_vk .* p[1:current_len]
end


function fft!(v; twiddles)
    n = length(v)
    @assert is_pow_2(n)

    fft_twiddles_parallel!(v; twiddles)
end

function ifft!(v; twiddles)
    n = length(v)
    @assert is_pow_2(n)

    ifft_twiddles!(v; twiddles)
end