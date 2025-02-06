@test function test_carryless_64(; n_iter = 100_000_000)
    a = rand(UInt64)
    t_begin = time()
    b = UInt128(1)
    for _ in 1:n_iter
        b ⊻= carryless_mul(UInt64(b&typemax(UInt64)), a)
    end
    t_end = time()
    @info "Throughput: $(n_iter/(t_end - t_begin))"
    return b
end


function test_carryless_128(; n_iter = 100_000_000)
    a = rand(UInt128)
    t_begin = time()
    for _ in 1:n_iter
        a ⊻= carryless_mul(a, a)[1]
    end
    t_end = time()
    @info "Throughput: $(n_iter/(t_end - t_begin))"
    return a
end

function test_mul_128(; n_iter = 100_000_000)
    a = GF2_128Elem(rand(UInt128))
    t_begin = time()
    for _ in 1:n_iter
        a += a*a
    end
    t_end = time()
    @info "Throughput: $(n_iter/(t_end - t_begin))"
end

function test_mat_mul(n)
    A = rand(GF2_128Elem, n, n)
    v = rand(GF2_128Elem, n)
    t_begin = time()
    y = A*v
    t_end = time()
    @info "Performed $(n^2) multiplies in $(t_end - t_begin) seconds"
    @info "Throughput (multiplies): $(n^2/(t_end - t_begin))"
    return y
end


function test()
    a = GF2_128Elem(1 << 24)
    b = GF2_128Elem(1)
    @show n_iter = 200_000
    for _ in 1:n_iter
        b *= a
    end
    return b
end


# XXX: Actually turn these into tests
