@testset "AAA" begin
    pts = 10 .^ range(-15,0,500); pts = [-reverse(pts); 0; pts]
    approx(f; kw...) = aaa(f; kw...)

    @testset "Basic functions" begin
        f = x -> sin(1/(1.05-x)); @test pass(f, approx(f), pts, atol=2e-13)
        f = x -> exp(-1/(x^2)); @test pass(f, approx(f), pts, rtol=4e-13)
        f = x -> exp(-100x^2); @test pass(f, approx(f), pts, rtol=2e-13)
        f = x -> exp(-10/(1.2-x)); @test pass(f, approx(f), pts, rtol=1e-12)
        f = x -> 1/(1+exp(100*(x+.5))); @test pass(f, approx(f), pts, atol=2e-13)
        f = x -> sin(100x) * exp(-10x^2); @test pass(f, approx(f), pts, atol=1e-11)
        f = x -> abs(x);  @test pass(f, approx(f), pts, atol=1e-8)
        f = x -> abs(x - 0.95);  @test pass(f, approx(f), pts, atol=1e-6)
        f = x -> cis(x); r = approx(f);
        @test sum(@. abs(r(pts))-1)/length(pts) < 2e-13
    end

    @testset "Fully discrete" begin
        f = x -> sin(1/(1.05-x)); @test pass(f, aaa(f.(pts), pts), pts, atol=2e-13)
        f = x -> exp(-1/(x^2)); @test pass(f, aaa(f.(pts), pts), pts, rtol=4e-13)
        f = x -> exp(-20x^2); @test pass(f, aaa(f.(pts), pts), pts, rtol=2e-13)
        f = x -> exp(-10/(1.2-x)); @test pass(f, aaa(f.(pts), pts), pts, rtol=1e-12)
        f = x -> 1/(1+exp(10*(x+.5))); @test pass(f, aaa(f.(pts), pts), pts, atol=2e-13)
        f = x -> sin(100x) * exp(-10x^2); @test pass(f, aaa(f.(pts), pts), pts, atol=1e-11)
    end

    @testset "Low-accuracy" begin
        f = x -> exp(3x);
        @test !pass(f, approx(f, tol=1e-4), pts, atol=1e-8)
    end

    @testset "Poles, zeros, residues" begin
        f = z -> (z+1) * (z+2) / ((z+3) * (z+4))
        r = approx(f)
        pol = poles(r)
        zer = roots(r)
        @test isapprox(sum(pol+zer), -10, atol=1e-12)

        f = x -> 2/(3 + x) + 5/(x - 2im);  r = approx(f)
        @test isapprox(prod(residues(r)[2]), 10, atol=1e-8)

        f = x -> sinpi(10x);  r = approx(f);
        @test isapprox(sort(abs.(roots(r)))[19], 0.9, atol=1e-12)

        f = z -> (z - (3 + 3im))/(z + 2);  r = approx(f)
        pol,zer = poles(r), roots(r)
        @test isapprox(pol[1]*zer[1], -6-6im, atol=1e-12)

        f = x -> exp(exp(x))/(x - 0.2im); r = approx(f);
        @test minimum(abs.(poles(r) .- .2im)) < 1e-12
    end

    @testset "Scaling" begin
        f = x -> 1e100*sin(x); @test pass(f, approx(f), pts, rtol=2e-13)
        f = x -> 1e-100*cos(x); @test pass(f, approx(f), pts, rtol=2e-13)
        x = range(-1, 1, 4000)
        r = aaa(sin.(10x), 1e100*x)
        @test r.(1e100x) ≈ sin.(10x) rtol=1e-13
        r = aaa(sin.(10x), 1e-100*x)
        @test r.(1e-100x) ≈ sin.(10x) rtol=1e-13
    end

    @testset "Specified " begin
        f = x -> 0; @test pass(f, approx(f, max_degree=0), pts, atol=2e-13)
        f = x -> x; @test pass(f, approx(f, max_degree=1), pts, atol=2e-13)
        f = x -> 1im*x; @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
        f = x -> x+x^2; @test pass(f, approx(f, max_degree=2), pts, atol=2e-13)
        f = x -> x+x^3; @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
        f = x -> 1/(1.1+x); @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
        f = x -> 1/(1+1im*x); @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
        f = x -> 1/(3+x+x^2); @test pass(f, approx(f, max_degree=2), pts, atol=2e-13)
        f = x -> 1/(1.01+x^3); @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
        f = x -> tanh(100x); @test pass(f, approx(f), pts, atol=2e-13)
        f = x -> tanh(100*(x-.2)); @test pass(f, approx(f), pts, atol=2e-13)
        f = x -> exp(x); @test pass(f, approx(f, tol=1e-13), pts, atol=2e-13)
        f = x -> cis(x); @test pass(f, approx(f), pts, atol=2e-13)
     end
end
