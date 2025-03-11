using SIMD

function mul(a::UInt64, b::UInt64)
    pmull_res = Vec(ccall("llvm.aarch64.neon.pmull64",
                          llvmcall,
                          NTuple{16, VecElement{UInt8}},
                          (UInt64,UInt64),
                          a, b))

    return reinterpret(UInt128, pmull_res)
end


function batch_16_poly_mul(a0::UInt16, a1::UInt16, b0::UInt16, b1::UInt16) 
    a1_sh = UInt64(a1)
    a1_sh <<= 32
    a = a0 ⊻ a1_sh 

    b1_sh = UInt64(b1)
    b1_sh <<= 32
    b = b0 ⊻ b1_sh 

    res = mul(a, b)
    a0b0 = convert(UInt32, res & typemax(UInt32))
    a1b1 = UInt32(res >> 64)

    @assert a0b0 == mul(UInt64(a0), UInt64(b0))
    @assert a1b1 == mul(UInt64(a1), UInt64(b1))

    return (a0b0, a1b1)
end

# irr(X) = X^16 + X^5 + X^3 + X^2 + 1
function simple_mul(a::UInt16, b::UInt16)
    res = mul(UInt64(a), UInt64(b))

    lo = convert(UInt16, res & typemax(UInt16))
    hi = UInt16(res >> 16)

    tmp = hi ⊻ (hi >> (16 - 5)) ⊻ (hi >> (16 - 3)) ⊻ (hi >> (16 - 2))
    res = lo ⊻ tmp ⊻ (tmp << 2) ⊻ (tmp << 3) ⊻ (tmp << 5)

    return res
end
