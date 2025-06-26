@testset "Arithmetic with $method" verbose=true for method in (Barycentric, Thiele)
    e = approximate(exp, unit_interval; method)
    t = approximate(tan, unit_interval; method)
    @test (e / e) ≈ 1
    @test (3im * t - 2im * t) ≈ 1im * t
    c = cos
    @testset "$(op)" for op in (+, -, *, /)
        @test values(op(e, 3.14im)) ≈ op.(values(e), 3.14im)
        @test nodes(op(e, 3.14im)) ≈ nodes(e)
        @test op(e, 3.14) ≈ z -> op(e(z), 3.14)
        @test op(3.14, e) ≈ z -> op(3.14, e(z))
        @test op(e, 3.14 + 2.72im) ≈ z -> op(e(z), 3.14 + 2.72im)
        @test op(3.14 + 2.72im, e) ≈ z -> op(3.14 + 2.72im, e(z))
        @test op(t, e) ≈ z -> op(t(z), e(z))
        @test op(c, e) ≈ z -> op(c(z), e(z))
        @test op(e, c) ≈ z -> op(e(z), c(z))
    end
    ec = approximate(exp, unit_circle; method)
    @test_throws DomainError ec + e
    @test !(ec ≈ e)
    cc = approximate(cis, unit_circle; method)
    @test ec + cc ≈ z -> exp(z) + cis(z)
end
