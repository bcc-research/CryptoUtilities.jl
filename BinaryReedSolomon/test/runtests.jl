using BinaryReedSolomon
using BinaryFields
using Test

@testset "FFT -> IFFT -> FFT" begin
    k = 10
    for T in [BinaryElem16, BinaryElem32, BinaryElem128]
        twiddles = compute_twiddles(T, k)
        v = rand(T, 2^k)
        v_fft = copy(v)
        fft!(v_fft; twiddles)
        @test v_fft != v
        ifft!(v_fft; twiddles)
        @test v_fft == v
    end
end

@testset "Reed-Solomon systematic" begin
end
