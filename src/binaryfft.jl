using Base.Threads

is_pow_2(n) = 2^(round(Int, log2(n))) == n

export mul_inplace!

"""
    fft_twiddles!(v; twiddles, idx=1)

Recursively apply twiddle factors to perform a Fast Fourier Transform (FFT) step in-place.

This function is a recursive implementation of a portion of a Fast Fourier Transform
algorithm, specifically designed to be used with pre-calculated twiddle factors.
It operates in-place on the input vector `v`, modifying it directly.

# Arguments
- `v`: A vector to be transformed in-place. Must have a length that is a power of 2.
- `twiddles`: A vector of pre-calculated twiddle factors. The structure and length
  of `twiddles` should be consistent with the intended FFT algorithm.
- `idx`: An index to select the appropriate twiddle factor from the `twiddles` vector
  for the current step. Defaults to 1 for the initial call.

# Keyword Arguments
- `twiddles`:  Required, the collection of twiddle factors.
- `idx`: Optional, the index into the twiddle factor collection, defaults to 1.

# Notes
- This function is designed to be used as a building block in a larger FFT implementation.
- The length of `v` is expected to be a power of 2 for the FFT algorithm to function correctly.
- The `idx` parameter is used to navigate the `twiddles` vector during the recursive calls,
  effectively selecting the correct twiddle factors for each stage of the FFT.
- The function modifies `v` in-place and does not return a new vector.
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

"""
    fft_twiddles_parallel!(v; twiddles, idx=1, thread_depth=nothing)

Recursively apply twiddle factors in parallel to perform a Fast Fourier Transform (FFT) step in-place.

This function is a parallel recursive implementation of a portion of a Fast Fourier Transform
algorithm, designed to be used with pre-calculated twiddle factors. It operates
in-place on the input vector `v`, modifying it directly, and leverages Julia's
threading capabilities for parallel execution.

# Arguments
- `v`: A vector to be transformed in-place. Must have a length that is a power of 2.
- `twiddles`: A vector of pre-calculated twiddle factors. The structure and length
  of `twiddles` should be consistent with the intended FFT algorithm.
- `idx`: An index to select the appropriate twiddle factor from the `twiddles` vector
  for the current step. Defaults to 1 for the initial call.
- `thread_depth`:  Controls the depth of parallel execution. If `nothing` (default), it is
  set to `round(Int, log2(nthreads()))` at the first call.  In subsequent recursive calls,
  it is decremented. When `thread_depth` reaches 0, the computation becomes serial,
  falling back to the sequential `fft_twiddles!` logic for the base cases.

# Keyword Arguments
- `twiddles`:  Required, the collection of twiddle factors.
- `idx`: Optional, the index into the twiddle factor collection, defaults to 1.
- `thread_depth`: Optional, the depth of threading recursion. If not provided, it's
  automatically determined based on the number of available threads.

# Notes
- This function is a parallel version of `fft_twiddles!` and is designed to exploit
  multi-core processors for faster FFT computations, especially for larger input sizes.
- The length of `v` is expected to be a power of 2 for the FFT algorithm to function correctly.
- The `idx` parameter is used to navigate the `twiddles` vector during the recursive calls,
  effectively selecting the correct twiddle factors for each stage of the FFT.
- The `thread_depth` parameter is crucial for controlling the granularity of parallelism.
  Setting it appropriately can help optimize performance by balancing threading overhead
  with the benefits of parallel computation. A depth that is too high might lead to excessive
  threading overhead, while a depth that is too low might not fully utilize available resources.
- The function modifies `v` in-place and does not return a new vector.
- If `thread_depth` is not explicitly provided in the initial call, the function will infer
  a default value based on the number of available threads and log this setting for informational
  purposes.
"""
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

"""
    split_half(v)

Splits a vector `v` into two halves.

Given a vector `v`, this function divides it into two equal parts and returns
them as two separate views.  If the length of `v` is `n`, the first view
will contain elements from index 1 to `n/2`, and the second view will contain
elements from index `n/2 + 1` to `n`.

# Arguments
- `v`: The input vector to be split.

# Returns
- A tuple containing two views of the input vector `v`. The first element
  of the tuple is a view of the first half of `v`, and the second element
  is a view of the second half of `v`.

# Notes
- This function returns views, not copies, of the original vector. This means
  that modifications to the returned views will directly affect the original
  vector `v`.
- It is assumed that the length of `v` is an even number, or that integer division
  (`div(n, 2)`) correctly represents the midpoint for the intended split.
- This function is often used in divide-and-conquer algorithms, such as the
  Fast Fourier Transform (FFT), where vectors are recursively split into halves.
"""
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
function mul_inplace!(v, λ)
    u, w = split_half(v)
    @views begin
        @. u += λ*w
        w .+= u
    end
end


next_s(s_prev, s_prev_at_root) = s_prev*s_prev + s_prev_at_root*s_prev

# s_i(v_i) = (s_{i-1}(v_i))^2 - s_{i-1}(v_{i - 1})*s_{i-1}(v_i)
# note that first two elemets of each layer are: 
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
    l0i += GF2_128Elem(bitreverse(UInt128(i)))
    layer[i] = l0i
  end

  # s0(v0)
  return GF2_128Elem(1)
end

function layer_i!(layer, layer_len, s_prev_at_root)
  prev_layer_len = 2 * layer_len
  s_at_root = compute_s_at_root(layer, s_prev_at_root)
  for (idx, s_prev) in enumerate(layer[1:2:prev_layer_len])
    layer[idx] = next_s(s_prev, s_prev_at_root)
  end
  return s_at_root
end 

function compute_twiddles(beta, k)
  twiddles = Vector{GF2_128Elem}(undef, 2^k - 1)
  layer = Vector{GF2_128Elem}(undef, 2^(k - 1))

  write_at = 2^(k - 1)
  s_prev_at_root = layer_0!(layer, beta, k)
  twiddles[write_at:end] .= layer 

  for _ in 1:(k - 1) 
    write_at = div(write_at, 2)
    # notice that layer_len = write_at 
    layer_len = write_at
    s_prev_at_root = layer_i!(layer, layer_len, s_prev_at_root)

    s_inv = inv(s_prev_at_root)
    twiddles[write_at:write_at+layer_len - 1] .= s_inv .* layer[0:layer_len]
  end

  return twiddles
end