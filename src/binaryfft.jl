using Base.Threads

export mul_inplace!

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

    mul_inplace!(v, twiddles[idx])
    
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

  ifft_matrix_inplace!(v, twiddles[idx])
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
          @info "Setting thread depth to  (did you launch julia with `--threads [thread_count]`?)"
        end
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

function split_half(v)
    n = length(v)
    n_div2 = div(n, 2)
    return @views v[1:n_div2], v[n_div2+1:end]
end

"""
    mul_inplace!(v, λ)

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
function mul_inplace!(v, λ)
    u, w = split_half(v)
    @views begin
        @. u += λ*w
        w .+= u
    end
end

# TODO: Add nice comments here

function ifft_matrix_inplace!(v, λ)
  lo, hi = split_half(v)
  @views begin
    hi .+= lo 
    @. l0 += λ*hi
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
    for i in 1:2^(k-1)
      l0i = beta
      l0i += GF2_128Elem((i-1) << 1)
      layer[i] = l0i
    end

    # s0(v0)
    return GF2_128Elem(1)
end

function layer_i!(layer, layer_len, s_prev_at_root)
    prev_layer_len = 2 * layer_len
    s_at_root = compute_s_at_root(layer, s_prev_at_root)
    for (idx, s_prev) in enumerate(@views layer[1:2:prev_layer_len])
      layer[idx] = next_s(s_prev, s_prev_at_root)
    end
    return s_at_root
end 

function compute_twiddles(beta, k)
    twiddles = Vector{GF2_128Elem}(undef, 2^k - 1) # 1 2 3 4 5 6 7
    layer = Vector{GF2_128Elem}(undef, 2^(k - 1))

    write_at = 2^(k - 1)
    s_prev_at_root = layer_0!(layer, beta, k)
    @views twiddles[write_at:end] .= layer 

    for _ in 1:(k - 1) 
      write_at >>= 1
      # notice that layer_len = write_at 
      layer_len = write_at
      s_prev_at_root = layer_i!(layer, layer_len, s_prev_at_root)

      s_inv = inv(s_prev_at_root)
      @views @. twiddles[write_at:write_at+layer_len-1] = s_inv * layer[1:layer_len]
    end

    return twiddles
end

is_pow_2(n) = 2^(round(Int, log2(n))) == n

export fft!

function fft!(v; twiddles=nothing, beta=nothing)
    n = length(v)
    @assert is_pow_2(n)
    k = round(Int, log2(n))

    beta = isnothing(beta) ? GF2_128Elem(0) : beta
    twiddles = isnothing(twiddles) ? compute_twiddles(beta, k) : twiddles

    fft_twiddles_parallel!(v; twiddles)
end

export ifft! 

function ifft!(v; twiddles=nothing, beta=nothing)
  n = length(v)
  @assert is_pow_2(n)
  k = round(Int, log2(n))

  beta = isnothing(beta) ? GF2_128Elem(0) : beta
  twiddles = isnothing(twiddles) ? compute_twiddles(beta, k) : twiddles

  ifft_twiddles_parallel!(v; twiddles)
end