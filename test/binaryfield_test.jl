@testset "BinaryField Addition" begin
    a = GF2_128Elem(rand(UInt128))
    zero_elem = zero(GF2_128Elem)

    @test a + a == zero_elem
    @test a + zero_elem == a
    @test zero_elem + a == a
end

@testset "BinaryField Multiplication" begin
    a = GF2_128Elem(rand(UInt128))
    one_elem = GF2_128Elem(UInt128(1))
    zero_elem = zero(GF2_128Elem)

    @test a * one_elem == a
    @test one_elem * a == a
    @test a * zero_elem == zero_elem
    @test zero_elem * a == zero_elem
end

@testset "BinaryField Inversion" begin
    one_elem = GF2_128Elem(1)
    zero_elem = zero(GF2_128Elem)

    # Test inverse for many random elements
    for _ in 1:100
        a = GF2_128Elem(rand(UInt128))
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