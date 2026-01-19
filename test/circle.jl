@testset "Circles" verbose=true begin
    pts = cispi.(2 * (0:500) / 500)
    UD = unit_disk
    UC = unit_circle

    @testset "Unit disk for $method" for method in (Barycentric, Thiele)
        approx(f; kw...) = approximate(f, UD; method, kw...)
        f = z -> sin(10z) * exp(-z^2); @test pass(f, approx(f), pts, rtol=2e-11)
        f = z -> sin(1/(1.1 - z)); @test pass(f, approx(f), pts, rtol=2e-13)
        f = sec; @test pass(f, approx(f, max_iter=15), pts, rtol=1e-6)
        f = z -> cos(sin(z)) + exp(z)/(z-1.1); @test pass(f, approx(f), pts, rtol=2e-13)
        f = x -> cis(x);  @test pass(f, approx(f), pts, atol=6e-13)
    end

    @testset "Unit circle for $method" for method in (Barycentric, Thiele)
        f = z -> abs(z-1im);  @test pass(f, approximate(f, UC), pts, rtol=2e-10)
        f = z -> tan(π*z);  @test pass(f, approximate(f, UC), pts, rtol=2e-13)
        f = z -> tanh(100z); @test pass(f, approximate(f, UC), pts, rtol=2e-13)
    end

    @testset "Float type conversion for $method" for method in (Barycentric, Thiele)
        f = z -> abs(z - 1im)
        r = approximate(f, UC; method)
        r32 = convert(Float32, r.fun)
        @test r32 isa method{Float32,ComplexF32}
    end

    @testset "Translate and scale for $method" for method in (Barycentric, Thiele)
        f = z -> sin(10z) * exp(-z^2)
        for (a, c) in ( (2.5, 0), (1, -1im), (0.4, -2))
            F = approximate(f, a*UC + c)
            @test pass(f, F, a*pts .+ c)
        end
        f = z -> 1e100sin(z); @test pass(f, approximate(f, UD), pts, rtol=2e-13)
        @test pass(f, approximate(f, UD, max_iter=12), pts, rtol=1e-6)
    end

    @testset "Poles, zeros, residues in $T for $method" for T in (Float64, Double64), method in (Barycentric, Thiele)
        UC = Circle{T}(0, 1)
        UD = interior(UC)
        f = z -> tan(T(π)*z);  F = approximate(f, UC)
        pol = poles(F); @test sort(abs.(pol))[1:5] ≈ 0.5*[1;1;3;3;5] atol=1e-3

        f = z -> exp(exp(z)) / (z - 1im // 5); pol = poles(approximate(f, UC));
        @test minimum(@. abs(pol - 1im // 5)) < 1000eps(T)

        f = z -> (z+1) * (z+2) / ((z+3) * (z+4));  F = approximate(f, UC)
        pol = poles(F);  zer = roots(F);
        @test isapprox(sum(pol+zer), -10, atol=1000eps(T))

        f = z -> 2/(3+z) + 5im / (z-2im);  F = approximate(f, UD)
        @test isapprox( prod(residues(F)[2]), 10im, atol=sqrt(eps(T)) )

        f = z -> (z-(3+3im))/(z+2);  F = approximate(f, UD)
        pol, zer = poles(F), roots(F);  @test isapprox(pol[1]*zer[1], -6-6im, atol=1000eps(T))
    end

    @testset "Tolerance" begin
        f = z -> exp(3*z); @test !pass(f, approximate(f, UD, tol=1e-4), pts, atol=1e-8)
        f = z -> exp(3*z); @test pass(f, approximate(f, UD, tol=1e-10), pts, atol=1e-8)
    end

    @testset "Low degree" begin
        f = x -> 0; @test pass(f, approximate(f, max_iter=1, UD), pts, atol=2e-13)
        f = x -> x; @test pass(f, approximate(f, max_iter=2, UD), pts, atol=2e-13)
        f = x -> x+x^2; @test pass(f, approximate(f, max_iter=3, UD), pts, atol=2e-13)
        f = x -> x+x^3; @test pass(f, approximate(f, max_iter=4, UD), pts, atol=2e-13)
        f = x -> x+x^3; @test !pass(f, approximate(f, max_iter=3, UD), pts, atol=2e-13)
        f = x -> 1/(3im + x + x^2); @test pass(f, approximate(f, max_iter=3, UC), pts, rtol=2e-13)
        f = x -> 1/(3im + x + x^2); @test !pass(f, approximate(f, max_iter=2, UC), pts, rtol=2e-13)
        f = x -> 1/(1.01 + x^3); @test pass(f, approximate(f, max_iter=4, UD), pts, rtol=2e-13)
        f = x -> 1/(1.01 + x^3); @test !pass(f, approximate(f, max_iter=3, UD), pts, rtol=2e-13)
    end
end
