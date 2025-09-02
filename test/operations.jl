@testset "Operations" begin
    @testset "Derivatives for $method" verbose=true for method in (Thiele, )
        using ComplexRegions
        @testset "Domain $it_d" for (it_d, domain) in enumerate((unit_interval, unit_disk, Shapes.square))
            @testset "Function $it_f" for (it_f, (f, df)) in enumerate((
                (x -> exp(x), x -> exp(x)),
                (x -> exp(-x), x -> -exp(-x)),
                (x -> cis(x), x -> 1im * cis(x)),
                (x -> x, x -> 1),
                (x -> 1im * x^2, x -> 2im * x),
                (x -> 1 / (1.1 - x), x -> 1 / (1.1 - x)^2),
                (x -> log(1.1 - x), x -> -1 / (1.1 - x)),
                (sin, cos),
            ))
                r = approximate(f, domain; method)
                @test isapprox(derivative(r; allowed=true), df)
                @test derivative(r.fun, 0.1) ≈ df(0.1)
            end
        end
    end

    @testset "Arithmetic with $method" verbose=true for method in (Barycentric, Thiele)
        @testset "Domain $iter" for (iter, domain) in enumerate((unit_interval, Shapes.square))
            e = approximate(exp, domain; method)
            t = approximate(tan, domain; method)
            @test (e / e) ≈ 1
            @test (3im * t - 2im * t) ≈ 1im * t
            c = cis
            ec = approximate(exp, unit_circle; method)
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
                @test_throws DomainError op(ec, e)
                @test !(ec ≈ e)
            end
        end
    end

    @testset "Arithmetic with zero for $method" verbose=true for method in (Barycentric, Thiele)
        r = approximate(exp, unit_interval; method)
        @test r + 0 ≈ r
        @test r - 0 ≈ r
        @test r * 0 ≈ 0
        @test_throws DomainError r / 0
        @test 0 + r ≈ r
        @test 0 - r ≈ -r
        @test 0 * r ≈ 0
        @test 0 / r ≈ 0
    end
end
