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