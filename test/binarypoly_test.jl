using CryptoUtilities
using Random
using Test

@testset "BinaryPoly" begin
    @testset "Conversion" begin
        @testset "Widening Conversions" begin
            poly8 = BinaryPoly8(0x0f)
            poly16 = convert(BinaryPoly16, poly8) # Widening: BinaryPoly8 -> BinaryPoly16
            @test poly16 isa BinaryPoly16
            @test CryptoUtilities.binary_val(poly16) === UInt16(0x0f)

            poly16 = BinaryPoly16(0x1234)
            poly32 = convert(BinaryPoly32, poly16) # Widening: BinaryPoly16 -> BinaryPoly32
            @test poly32 isa BinaryPoly32
            @test CryptoUtilities.binary_val(poly32) === UInt32(0x1234)

            poly64 = BinaryPoly64(0x1234567890abcdef)
            poly128 = convert(BinaryPoly128, poly64) # Widening: BinaryPoly64 -> BinaryPoly128
            @test poly128 isa BinaryPoly128
            @test CryptoUtilities.binary_val(poly128) === UInt128(0x1234567890abcdef)
        end

        @testset "Narrowing Conversions - Exact Fit" begin
            poly128_small = BinaryPoly128(UInt128(0x0f))
            poly8_from_small_128 = convert(BinaryPoly8, poly128_small) # Narrowing, value fits: BinaryPoly128 -> BinaryPoly8
            @test poly8_from_small_128 isa BinaryPoly8
            @test CryptoUtilities.binary_val(poly8_from_small_128) === UInt8(0x0f)

            poly64_ff = BinaryPoly64(UInt64(0xff))
            poly8_from_64_ff = convert(BinaryPoly8, poly64_ff) # Narrowing, value fits in UInt8's range
            @test poly8_from_64_ff isa BinaryPoly8
            @test CryptoUtilities.binary_val(poly8_from_64_ff) == UInt8(0xff)
        end

        @testset "Narrowing Conversions - Error on potential data loss" begin
            poly128_large = BinaryPoly128(UInt128(0x1234567890abcdef0000000000000000))
            @test_throws InexactError convert(BinaryPoly64, poly128_large) # Narrowing, potential data loss: BinaryPoly128 -> BinaryPoly64

            poly64_100 = BinaryPoly64(UInt64(0x100)) # Value too large for UInt8 (256)
            @test_throws InexactError convert(BinaryPoly8, poly64_100) # Narrowing, potential data loss: BinaryPoly64 -> BinaryPoly8
        end
    end

    @testset "Addition" begin
        for BinaryPolyType in [BinaryPoly8, BinaryPoly16, BinaryPoly32, BinaryPoly64, BinaryPoly128]
            UIntType = CryptoUtilities.primitive_type(BinaryPolyType)
            val1 = rand(UIntType)
            val2 = rand(UIntType)
            poly1 = BinaryPolyType(val1)
            poly2 = BinaryPolyType(val2)

            sum_poly = poly1 + poly2
            @test sum_poly isa BinaryPolyType
            @test CryptoUtilities.binary_val(sum_poly) === val1 ⊻ val2

            @test zero(BinaryPolyType) + poly1 === poly1
            @test poly1 + zero(BinaryPolyType) === poly1
            @test poly1 + poly1 === zero(BinaryPolyType) # a + a = 0 in GF(2)
        end
    end

    @testset "Multiplication BinaryPoly64" begin
        poly_a = BinaryPoly64(UInt64(0b10110010_01010110_11001111_00010010_10010011_01101010_00001111_11110000))
        poly_b = BinaryPoly64(UInt64(0b01101001_10100011_00110101_10001100_01011010_00110011_11110000_00001111))

        prod_poly = poly_a * poly_b
        @test prod_poly isa BinaryPoly128
        # Value from sage
        @test CryptoUtilities.binary_val(prod_poly) === UInt128(0b111111101100010010001010011110101000011010001101101110110100101011010010001011010000101010000101110111010100110101000001010000)
    end
end
