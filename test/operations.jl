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
                r = approximate(f, domain, method())
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
            e = approximate(exp, domain, method())
            t = approximate(tan, domain, method())
            @test (e / e) ≈ 1
            @test (3im * t - 2im * t) ≈ 1im * t
            c = cis
            ec = approximate(exp, unit_circle, method())
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
        r = approximate(exp, unit_interval, method())
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

@testset "QTCF finite-block construction" begin
    @testset "Rejection before mutation with $method" for method in (Classic(), OneDiv())
        g = QTCF([0.0], [1.0]; method)
        before = copy(g)
        for (z, y) in (([1.0], [1.0]), ([1.0, 2.0], [1.0, 1.0]),
                       ([1.0, 2.0], [1.0, 2.0]), ([1.0], [Inf]), ([1.0], [NaN]))
            @test_throws ArgumentError RFA.add_node!(g, z, y, method)
            for name in fieldnames(typeof(g))
                @test isequal(getfield(g, name), getfield(before, name))
            end
        end
        @test_throws ArgumentError QTCF([0.0, 1.0], [1.0, 1.0]; method)
        @test_throws ArgumentError QTCF([0.0, 1.0, 2.0], [1.0, 1.0, 1.0];
            block_sizes=[1, 2], method)
        @test_throws ArgumentError QTCF([0.0, 1.0], [1.0, 1.0], [[1.0], [Inf]])
        @test_throws RFA.NaNException QTCF([0.0], [1.0], [[NaN]])
        @test_throws ArgumentError QTCF([0.0], [Inf])
        @test_throws ArgumentError QTCF{Float32}([0.0], [1e100])
        @test_throws ArgumentError QTCF{Float32}([0.0], [1.0], [[1e100]])
        @test_throws ArgumentError QTCF{Float32}([1.0, 1.0 + 1e-8], [1.0, 2.0])
        @test RFA.add_node!(g, 1.0, 2.0, method) === g
        @test all(isfinite, Iterators.flatten(weights(g)))
        @test evaluate(g, 0.5, OneDiv()) ≈ 1.5
    end

    @testset "Intermediate infinity and finite zero weights" begin
        for method in (Classic(), OneDiv())
            g = QTCF([0.0, 1.0], [1.0, 2.0]; method)
            @test isinf(RFA._qtcf_ratio(2.0, 1.0 - first(values(g))))
            @test RFA.add_node!(g, 2.0, 1.0, method) === g
            @test iszero(only(last(weights(g))))
            @test all(isfinite, Iterators.flatten(weights(g)))
            @test evaluate(g, 0.5, OneDiv()) == 1.0
        end
    end

    @testset "No previous-node attainment audit" begin
        for method in (Classic(), OneDiv())
            g = QTCF([0.0], [0.0]; method)
            @test RFA.add_node!(g, [-1.0, 1.0], [1.0, 1.0], method) === g
            @test g.P == g.Q == [0.0, 1.0]
            @test isnan(evaluate(g, 0.0, OneDiv()))
            @test evaluate(g, 0.5, Classic()) == 1.0
            @test RFA.add_node!(g, 0.5, 0.5, method) === g
        end
    end

    @testset "Finite coefficients in all approximation paths" begin
        z = collect(range(-1.0, 1.0, length=25))
        for reflected in (false, true), continuous in (false, true)
            kw = reflected ? (; reflection=x -> -x) : (;)
            r = continuous ? approximate(exp, Segment(-1.0, 1.0), QTCF();
                tol=1e-10, max_iter=8, stagnation=20, kw...) :
                approximate(exp.(z), z, QTCF(); tol=1e-10, max_iter=8,
                    stagnation=20, kw...)
            @test all(all(isfinite, Iterators.flatten(weights(h.interpolant))) for h in r.history)
            @test maximum(abs.(r.(z) .- exp.(z))) < 1e-6
            @test get_history(r; get_poles=false)[end] == status(r).best
            @test nodes(rewind(r, 1)) == nodes(first(r.history).interpolant)
        end
    end

    @testset "Distinct nodes under mixed residual selection" begin
        z = collect(range(-1.0, 1.0, length=9))
        kw = (; tol=0.0, max_iter=3, stagnation=20,
            first_residual=:true_residual, second_residual=:linearized_residual)
        fits = (approximate(exp.(z), z, QTCF(); kw..., noise_floor_factor=0),
                approximate(exp, Segment(-1.0, 1.0), QTCF(); kw...,
                    refinement=3, initial_refinement=3),
                approximate(z -> real(z), ClosedCurve(t -> cispi(4t)), QTCF();
                    tol=0.0, max_iter=3, stagnation=20, refinement=7,
                    initial_refinement=7))
        for r in fits
            g = last(r.history).interpolant
            @test length(g.pairs) > 1
            @test all(!RFA._same(a, b) for (a, b) in g.pairs[2:end])
        end
    end

    @testset "Rejected candidates do not prevent alternative blocks" begin
        r = approximate([0.0, 2.0, 0.0, 1.0, 1.0], [0.0, 1.0, -1.0, 2.0, -2.0],
            QTCF(); reflection=z -> -z, tol=0.0, noise_floor_factor=0,
            max_iter=2, stagnation=20)
        g = last(r.history).interpolant
        @test g.pairs == [(0.0, 0.0), (2.0, -2.0)]
        @test all(isfinite, Iterators.flatten(weights(g)))
        r = approximate([0.0, 0.0, 1.0, 2.0], [0.0, 1.0, 2.0, 3.0],
            QTCF(); tol=0.0, noise_floor_factor=0, max_iter=2, stagnation=20)
        @test last(r.history).interpolant.pairs == [(0.0, 0.0), (3.0, 2.0)]
    end

    @testset "Stable evaluation and derivative recovery" begin
        g = QTCF([0.0, 1.0, 2.0, 3.0], zeros(4), [1.0, 1e150, 1e150, 1e150])
        @test evaluate(g, 0.7, OneDiv()) ≈ 1.0
        v = derivative(g, [0, 1])(0.7)
        @test all(isfinite, v)
        @test v[1] ≈ 1.0
        @test isapprox(v[2], 1e-150; rtol=2e-14, atol=0)
        for T in (Float32, Float64, BigFloat)
            r = QTCF(T[-1, 0, 1], exp.(T[-1, 0, 1]))
            z = T[0.2, 0.4, 0.7]
            @test evaluate(r, z, OneDiv()) ≈ evaluate(r, z, Classic())
            @test evaluate(r, reshape(z, 1, 3), OneDiv()) ≈ reshape(r.(z), 1, 3)
        end
    end

    @testset "Contour residues at an exact repeated pole" begin
        for c in (0.0, 3.0)
            g = QTCF([1.0, c], [0.0, 0.0, 1.0], Tuple{Float64,Float64}[],
                Vector{Float64}[], Vector{Float64}[], Float64[], Float64[])
            p, res = residues(g)
            @test length(p) == 2
            @test all(iszero, p)
            @test all(isfinite, res)
            @test all(isapprox(v, c; atol=1e-12) for v in res)
            @test_throws ArgumentError convert(PartialFractions, g)
        end
    end
end
