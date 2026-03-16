all_types = [BinaryElem16, BinaryElem32, BinaryElem128]

@testset "BinaryField Addition" begin
    for T in all_types
        a = rand(T)
        zero_elem = zero(T)

        @test a + a == zero_elem
        @test a + zero_elem == a
        @test zero_elem + a == a
    end
end

@testset "BinaryField Multiplication" begin
    for T in all_types
        a = rand(T)
        zero_elem = zero(T)
        one_elem = one(T)

        @test a * one_elem == a
        @test one_elem * a == a
        @test a * zero_elem == zero_elem
        @test zero_elem * a == zero_elem
    end
end

@testset "BinaryField Inversion" begin
    for T in all_types
        one_elem = one(T)
        zero_elem = zero(T)

        # Test inverse for many random elements
        for _ in 1:100
            a = rand(T)
            if a != zero_elem
                @test a * inv(a) == one_elem
            end
        end

        # Test inverse of 1
        @test inv(one_elem) == one_elem

        # Inverse of zero is not defined, but the code might return something.
        # We should check if it throws an error or returns a specific value (though it's not mathematically defined).
        # For now, let's just check that it doesn't crash.
        @test_throws AssertionError inv(zero_elem) # div_irreducible asserts a != 0
    end
end

@testset "BinaryField Powers" begin
    for T in [BinaryElem16, BinaryElem32]
        for _ in 1:100
            a = rand(T)

            # Test base cases
            @test a^0 == one(T)
            @test a^1 == a

            # Test subfield condition
            @test a^(2^BinaryFields.bitsize(T)) == a
            if a != zero(T)
                @test a^(2^BinaryFields.bitsize(T) - 2)*a == one(T)
            end
        end
    end
end

@testset "BinaryField Convert and Mapping" begin
    @testset "promote_rule" begin
        for T in [BinaryElem16, BinaryElem32]
            for _ in 1:100
                a = rand(T)
                a_conv = convert(BinaryElem128, a)
                
                @test a_conv == a * BinaryElem128(1)
            end
        end
    end

    @testset "Subfield condition" begin
        for T in [BinaryElem16, BinaryElem32]
            for _ in 1:100
                a = rand(T)
                a_conv = convert(BinaryElem128, a)

                @test a_conv^(2^BinaryFields.bitsize(T)) == a_conv
                if a != zero(T)
                    @test a_conv^(2^BinaryFields.bitsize(T) - 2)*a_conv == one(BinaryElem128)
                end
            end
        end
    end

    @testset "Addition and multiplication preserved" begin
        for T in [BinaryElem16, BinaryElem32]
            for _ in 1:100
                a = rand(T)
                a_conv = convert(BinaryElem128, a)

                b = rand(T)
                b_conv = convert(BinaryElem128, b)

                @test convert(BinaryElem128, a+b) == a_conv+b_conv
                @test convert(BinaryElem128, a*b) == a_conv*b_conv
            end
        end
    end

    @testset "Inverse preserved" begin
        for T in [BinaryElem16, BinaryElem32]
            for _ in 1:100
                a = rand(T)
                a_conv = convert(BinaryElem128, a)

                if a != zero(T)
                    @test inv(a_conv) == convert(BinaryElem128, inv(a))
                end
            end
        end
    end

end

@testset "BinaryPoly divrem" begin
    for T in [BinaryPoly16, BinaryPoly32]
        for _ in 1:100
            a = rand(T)
            b = rand(T)
            if binary_val(b) == 0
                b = one(T)
            end
            q, r = BinaryFields.divrem(a, b)
            @test q * b + r == a
            @test BinaryFields.degree(r) < BinaryFields.degree(b) || binary_val(r) == 0
        end
    end
end
