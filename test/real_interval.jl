@testset "Real intervals" verbose=true begin
    domain = Dict()
    test_points = Dict()
    for T in (Float64, Double64)
        domain[T] = Segment{T}(-1, 1)
        z = T(10) .^ range(T(-15), T(0), 500);
        test_points[T] = [-reverse(z); 0; z]
    end

    @testset "Basic functions for $method" for method in (Barycentric, Thiele)
        T = Float64
        tol = 2000*eps(T)
        pts = test_points[T]
        approx(f; kw...) = approximate(f, domain[T]; method, kw...)
        @testset "Function $iter" for (iter, f) in enumerate((
            exp,
            cis,
            x -> 1im*x + exp(-1 / x^2),
            x -> 1 / (1.1 - x),
            x -> log(1.1 - x),
            x -> sin(1 / (21//20 - x)),
            x -> abs(x + 1//2 + 1im//100),
            x -> exp(-100x^2),
            x -> sin(80x) * exp(-10x^2),
            x -> 1im + exp(-10 / (6//5 - x)),
            x -> 10im*x + tanh(100*(x - 1//5)),
            x -> x + tanh(100x),
            ))
            r = approx(f)
            @test pass(f, r, pts, rtol=tol)
        end
        # alternate checking methods
        r = approx(cis)
        @test isapprox(r, cis; rtol=tol)
        _, err = check(r; quiet=true); @test maximum(abs, err) < tol
        @test sum(@. abs(r(pts))-1)/length(pts) < tol

    end

    @testset "Double64 for $method" for method in (Barycentric, Thiele)
        T = Double64
        tol = 2000*eps(T)
        pts = test_points[T]
        approx(f; kw...) = approximate(f, domain[T]; method, kw...)
        @testset "Function $iter" for (iter, f) in enumerate((
            x -> cis(x),
            x -> exp(x),
            x -> abs(x + 1//4 + 1im//100),
            x -> exp(-1 / x^2),
            x -> exp(-60*(x + 1//6)^2),
            x -> 1im + sin(1 / (21//20 - x)),
            x -> exp(-10 / (6//5 - x)),
            x -> tanh(100*(x - 1//5))
            ))
            r = approx(f)
            @test pass(f, r, pts, rtol=tol)
        end
    end

    @testset "Tolerance" begin
        T = Float64
        method = Barycentric
        pts = test_points[T]
        approx(f; kw...) = approximate(f, domain[T]; method, kw...)
        f = x -> exp(3x);
        r = approx(f, tol=1e-5)
        @test !pass(f, r, pts, atol=1e-7)
        @test pass(f, r, pts, atol=5e-5)
        f = x -> abs(x);  @test pass(f, approx(f, stagnation=30), pts, atol=1e-10)
        f = x -> abs(x - 0.95);  @test pass(f, approx(f, stagnation=30), pts, atol=1e-9)
    end

    @testset "Nodes, values, degree for Barycentric" begin
        r = approximate(exp, unit_interval; method=Barycentric)
        @test length(nodes(r)) == 6
        @test length(weights(r)) == 6
        @test minimum(nodes(r)) ≈ -1
        @test length(values(r)) == 6
        @test minimum(values(r)) ≈ exp(-1)
        @test degree(r) == length(nodes(r)) - 1
        @test degrees(r) == (degree(r), degree(r))
        @test isapprox(r, exp)
        deg, err, zp, allowed, best = get_history(r)
        @test deg[end] == degree(r)
        @test length(deg) == length(err) == length(allowed)
        r = approximate(exp, unit_interval; method=Barycentric, allowed=true)
        deg, err, zp, allowed, best = get_history(r)
        @test deg[end] == degree(r)
        @test length(deg) == length(err) == length(allowed)
    end

    @testset "Nodes, values, degree for Thiele" begin
        r = approximate(exp, unit_interval; method=Thiele)
        @test length(nodes(r)) == 11
        @test length(weights(r)) == 11
        @test minimum(nodes(r)) ≈ -1
        @test length(values(r)) == 11
        @test minimum(values(r)) ≈ exp(-1)
        @test degree(r) == 5
        @test degrees(r) == (5, 5)
        @test isapprox(r, exp)
        deg, err, zp, allowed, best = get_history(r)
        @test deg[end] == degree(r)
        @test length(deg) == length(err) == length(allowed)
        r = approximate(exp, unit_interval; method=Barycentric, allowed=true)
        deg, err, zp, allowed, best = get_history(r)
        @test deg[end] == degree(r)
        @test length(deg) == length(err) == length(allowed)
    end

    @testset "Poles, zeros, residues in $T for Barycentric" for T in (Float64, Double64)
        approx(f; kw...) = approximate(f, domain[T]; method=Barycentric, kw...)
        f = z -> (z+1) * (z+2) / ((z+3) * (z+4))
        r = approx(f)
        pol = poles(r)
        zer = roots(r)
        @test isapprox(sum(pol+zer), -10, rtol=5000*eps(T))

        f = x -> 2/(3 + x) + 5/(x - 2im);  r = approx(f)
        @test isapprox(prod(residues(r)[2]), 10, rtol=sqrt(eps(T)))

        f = x -> sinpi(10x);  r = approx(f);
        zer = roots(r)
        miss = minimum(abs, zer .- 9//10)
        @test miss < 1e3*eps(T)

        f = z -> (z - (3 + 3im))/(z + 2);  r = approx(f)
        pol,zer = poles(r), roots(r)
        @test isapprox(pol[1]*zer[1], -6-6im, rtol=5000eps(T))

        f = x -> exp(exp(x))/(x - 1im // 5); r = approx(f);
        @test minimum(abs.(poles(r) .- 1im // 5)) < 1000eps(T)
    end

    @testset "Poles, zeros, residues for Thiele" begin
        approx(f; kw...) = approximate(f, domain[Float64]; method=Thiele, kw...)
        f = z -> (z+1) * (z+2) / ((z+3) * (z+4))
        r = approx(f)
        pol = poles(r)
        zer = roots(r)
        @test isapprox(sum(pol + zer), -10, rtol=5000*eps())

        f = x -> 2im / (1.2 + x) + 3 / (x - 1.25im);  r = approx(f)
        @test isapprox(prod(residues(r)[2]), 6im, rtol=1e-13)

        f = x -> sinpi(10x) + x - 9//10;  r = approx(f);
        zer = roots(r)
        miss = minimum(abs, zer .- 9//10)
        @test miss < 1e-7

        f = z -> (z - (3 + 3im))/(z + 2);  r = approx(f)
        pol, zer = poles(r), roots(r)
        @test isapprox(pol[1]*zer[1], -6-6im, rtol=5000*eps())

    end

    @testset "Thiele direct construction" begin
        f = x -> sinpi(x) + x - 9//10;
        x = range(-1, 1, 20);
        r = Thiele(x, f.(x))
        @test maximum(abs(r(x) - f(x)) for x in range(-1,1,1000)) < 1e-10
    end

    @testset "Vertical scaling in $T" for T in (Float64, Double64)
        pts = test_points[T]
        approx(f; kw...) = approximate(f, domain[T]; method=Barycentric, kw...)
        f = x -> T(10)^50 * sin(x); @test pass(f, approx(f), pts, rtol=2000*eps(T))
        f = x -> T(10)^(-50) * cos(x); @test pass(f, approx(f), pts, rtol=2000*eps(T))
        # TODO: Horizontal scaling is broken, because pole computation uses 1
        #s = T(10)^50
        #f = x -> sin(s * x)
        #@test pass(f, approximate(f, domain[T] / s), pts / s, rtol=2000*eps(T))
    end

     @testset "Low degree" begin
        T = Float64
        tol = 2000*eps(T)
        pts = test_points[T]
        approx(f; kw...) = approximate(f, domain[T]; method=Barycentric, kw...)
        f = x -> 0; @test pass(f, approx(f, max_iter=1), pts, atol=tol)
        f = x -> x; @test pass(f, approx(f, max_iter=2), pts, atol=tol)
        f = x -> 1im*x; @test pass(f, approx(f, max_iter=2), pts, atol=tol)
        f = x -> x+x^2; @test pass(f, approx(f, max_iter=3), pts, atol=tol)
        f = x -> x+x^3; @test pass(f, approx(f, max_iter=4), pts, atol=tol)
        f = x -> 1/(1.1+x); @test pass(f, approx(f, max_iter=3), pts, atol=tol)
        f = x -> 1/(1+1im*x); @test pass(f, approx(f, max_iter=3), pts, atol=tol)
        f = x -> 1/(3+x+x^2); @test pass(f, approx(f, max_iter=3), pts, atol=tol)
        f = x -> 1/(1.01+x^3); @test pass(f, approx(f, max_iter=4), pts, atol=tol)
    end

    @testset "Interval [$a, $b]" for (a, b) in ((-2, 3), (0, 4), (-2e-4, 0), (-3e3, 5e6))
        pts = range(a, b, 1000)
        approx(f; kw...) = approximate(f, Segment(a, b); method=Barycentric, kw...)
        toler(f) = 1000 * eps() * max(b - a, abs(f(a)), abs(f(b)))
        f = x -> 1 / sin(a + (b-a)*1.05im - x); @test pass(f, approx(f), pts, atol=toler(f))
        f = x -> exp(-10/(a + 1.1*(b-a) - x)); @test pass(f, approx(f), pts, atol=toler(f))
    end
end
