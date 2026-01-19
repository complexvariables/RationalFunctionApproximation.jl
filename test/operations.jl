@testset "Operations" begin
    @testset "Derivatives for $method" verbose=true for method in (Thiele, Barycentric)
        using ComplexRegions
        @testset "Domain $it_d" for (it_d, domain) in enumerate((unit_interval, unit_disk, Shapes.square))
            @testset "Function $it_f" for (it_f, (f, df, d2f)) in enumerate((
                (exp, exp, exp),
                (x -> 3., x -> 0., x -> 0.),
                (x -> x, x -> 1, x -> 0),
                (x -> exp(-x), x -> -exp(-x), x -> exp(-x)),
                (x -> cis(x), x -> 1im * cis(x), x -> -cis(x)),
                (x -> 1im * x^2, x -> 2im * x, x -> 2im),
                (x -> 1 / (1.1 - x), x -> 1 / (1.1 - x)^2, x -> 2 / (1.1 - x)^3),
                (x -> log(1.1 - x), x -> -1 / (1.1 - x), x -> -1 / (1.1 - x)^2),
                (sin, cos, x -> -sin(x)),
            ))
                r = approximate(f, domain; method)
                @test isapprox(derivative(r; allowed=true), df, atol=sqrt(eps()))
                @test isapprox(derivative(r, 2; allowed=true), d2f, atol=50sqrt(eps()))
                vals = derivative(r.fun, 0:2)(0.25)
                @test isapprox(vals[1], f(0.25), atol=sqrt(eps()))
                @test isapprox(vals[2], df(0.25), atol=sqrt(eps()))
                @test isapprox(vals[3], d2f(0.25), atol=50sqrt(eps()))
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
                @test isapprox(op(c, e), z -> op(c(z), e(z)), atol=sqrt(eps()))
                @test isapprox(op(e, c), z -> op(e(z), c(z)), atol=sqrt(eps()))
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
