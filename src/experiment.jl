using SIMD

function mul(a::UInt64, b::UInt64)
    pmull_res = Vec(ccall("llvm.aarch64.neon.pmull64",
                          llvmcall,
                          NTuple{16, VecElement{UInt8}},
                          (UInt64,UInt64),
                          a, b))

    return reinterpret(UInt128, pmull_res)
end


function batch_16_poly_mul(a, b) 
    a0, a1 = a
    b0, b1 = b
    a1_sh = UInt64(a1)
    a1_sh <<= 32
    a = a0 ⊻ a1_sh 

    b1_sh = UInt64(b1)
    b1_sh <<= 32
    b = b0 ⊻ b1_sh 

    res = mul(a, b)
    a0b0 = convert(UInt32, res & typemax(UInt32))
    a1b1 = UInt32(res >> 64)

    # @assert a0b0 == mul(UInt64(a0), UInt64(b0))
    # @assert a1b1 == mul(UInt64(a1), UInt64(b1))

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

function other_mul(a)
    for _ in 1:1000
        a += a*a
    end
end

function batch_16_poly_mul_ll_simd(a0::UInt32, a1::UInt32, b0::UInt32, b1::UInt32)
    llvm_ir = """
    define { i64, i64 } @batch_mul(i32 %a0, i32 %a1, i32 %b0, i32 %b1) {
        ; Extend inputs to 64-bit and pack
        %a0_zext = zext i32 %a0 to i64
        %a1_zext = zext i32 %a1 to i64
        %a1_shl = shl i64 %a1_zext, 32
        %a = xor i64 %a0_zext, %a1_shl

        %b0_zext = zext i32 %b0 to i64
        %b1_zext = zext i32 %b1 to i64
        %b1_shl = shl i64 %b1_zext, 32
        %b = xor i64 %b0_zext, %b1_shl

        ; Create return struct: { i64, i64 }
        %result = insertvalue { i64, i64 } undef, i64 %a, 0
        %result2 = insertvalue { i64, i64 } %result, i64 %b, 1
        ret { i64, i64 } %result2
    }
    """
    return Base.llvmcall((llvm_ir, "batch_mul"), Tuple{UInt32,UInt32},
                         Tuple{UInt32,UInt32,UInt32,UInt32}, a0, a1, b0, b1)
end
