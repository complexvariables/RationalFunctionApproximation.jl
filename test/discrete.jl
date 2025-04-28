@testset "Discrete interval for $method" for method in (Barycentric, Thiele)
    T = Float64
    tol = 2000*eps(T)
    z = T(10) .^ range(T(-15), T(0), 500);
    pts = [-reverse(z); 0; z]
    approx(f; kw...) = approximate(f, pts; method, kw...)
    f = x -> abs(x + 1//2 + 1im//100); @test pass(f, approx(f; stagnation=30), pts; rtol=tol)
    f = x -> sin(1 / (21//20 - x)); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> 1im*x + exp(-1 / x^2); @test pass(f, approx(f; stagnation=30), pts; rtol=tol)
    if method != Thiele
        f = x -> exp(100x^2); @test pass(f, approx(f), pts; rtol=tol)
        f = x -> x + sin(80x) * exp(-10x^2); @test pass(f, approx(f; stagnation=30), pts; rtol=tol)
    end
    f = x -> exp(-10 / (6//5 - x)); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> 10im*x + tanh(100*(x - 1//5)); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> x + tanh(100x); @test pass(f, approx(f; stagnation=20), pts, rtol=tol)
    f = x -> exp(x); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> cis(x); @test pass(f, approx(f), pts, rtol=tol)
end


@testset "Discrete interval for Double64" begin
    method = Barycentric
    T = Double64
    tol = 2000*eps(T)
    z = T(10) .^ range(T(-15), T(0), 500);
    pts = [-reverse(z); 0; z]
    approx(f; kw...) = approximate(f, pts; method, kw...)
    f = x -> sin(1 / (21//20 - x)); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> 1im*x + exp(-1 / x^2); @test pass(f, approx(f; stagnation=30), pts; rtol=tol)
    f = x -> x + sin(80x) * exp(-10x^2); @test pass(f, approx(f; stagnation=30), pts; rtol=tol)
    f = x -> 10im*x + tanh(100*(x - 1//5)); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> exp(x); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> cis(x); @test pass(f, approx(f), pts, rtol=tol)
end

@testset "Discrete circle for $method" for method in (Barycentric, Thiele)
    pts = cispi.(2 * (0:1000) / 1000)
    approx(f; kw...) = approximate(f, pts; method, kw...)
    f = z -> sin(10z) * exp(-z^2); @test pass(f, approx(f), pts, rtol=2e-11)
    f = z -> sin(1/(1.1 - z)); @test pass(f, approx(f), pts, rtol=2e-13)
    f = sec; @test pass(f, approx(f, max_iter=15), pts, rtol=1e-6)
    f = z -> cos(sin(z)) + exp(z)/(z-1.1); @test pass(f, approx(f), pts, rtol=2e-13)
    f = x -> cis(x);  @test pass(f, approx(f), pts, atol=6e-13)
    f = z -> tan(π*z);  @test pass(f, approximate(f, UC), pts, rtol=2e-13)
    f = z -> tanh(100z); @test pass(f, approximate(f, UC), pts, rtol=2e-13)
end

@testset "Discrete interval, low accuracy" begin
    pts = range(-1, 1, 1001)
    approx(f; kw...) = approximate(f, pts; method=Barycentric, kw...)
    f = x -> exp(3x);
    r = approx(f, tol=1e-5)
    @test !pass(f, r, pts, atol=1e-9)
    @test pass(f, r, pts, atol=2e-5)
    f = x -> abs(x);  @test pass(f, approx(f), pts, atol=1e-7)
end

@testset "Poles, zeros, residues in $T" for T in (Float64,Double64)
    approx(f; kw...) = approximate(f, range(T(-1), T(1), 1001); method=Barycentric, kw...)
    f = z -> (z+1) * (z+2) / ((z+3) * (z+4))
    r = approx(f)
    pol = poles(r)
    zer = roots(r)
    @test isapprox(sum(pol+zer), -10, rtol=5000*eps(T))

    f = x -> 2/(3 + x) + 5/(x - 2im);  r = approx(f)
    @test isapprox(prod(values(residues(r))), 10, rtol=sqrt(eps(T)))

    f = x -> sinpi(10x);  r = approx(f);
    zer = roots(r)
    miss = minimum(abs, zer .- 9//10)
    @test miss < 1e3*eps(T)

    f = z -> (z - (3 + 3im))/(z + 2);  r = approx(f)
    pol,zer = poles(r), roots(r)
    @test isapprox(pol[1]*zer[1], -6-6im, rtol=5000*eps(T))
end

@testset "Vertical scaling in $T" for T in (Float64, Double64)
    pts = range(T(-1), T(1), 1001)
    approx(f; kw...) = approximate(f, pts; method=Barycentric, kw...)
    f = x -> T(10)^50 * sin(x); @test pass(f, approx(f), pts, rtol=2000*eps(T))
    f = x -> T(10)^(-50) * cos(x); @test pass(f, approx(f), pts, rtol=2000*eps(T))
end


@testset "Polynomials and reciprocals" begin
    pts = range(-1, 1, 1001)
    tol = 2000*eps(T)
    approx(f; kw...) = approximate(f, pts; method=Barycentric, kw...)
    f = x -> 0; @test pass(f, approx(f), pts, atol=tol)
    f = x -> x; @test pass(f, approx(f), pts, atol=tol)
    f = x -> 1im*x; @test pass(f, approx(f), pts, atol=tol)
    f = x -> x + x^3; @test pass(f, approx(f), pts, atol=tol)
    f = x -> x + x^2; @test pass(f, approx(f), pts, atol=tol)
    f = x -> 1/(1.1 + x); @test pass(f, approx(f), pts, atol=tol)
    f = x -> 1/(1 + 1im*x); @test pass(f, approx(f), pts, atol=tol)
    f = x -> 1/(3 + x + x^2); @test pass(f, approx(f), pts, atol=tol)
    f = x -> 1/(1.01 + x^3); @test pass(f, approx(f), pts, atol=tol)
end

@testset "Specified degree" begin
    tol = 2000*eps(T)
    pts = range(-1, 1, 1001)
    approx(f; kw...) = approximate(f, pts; method=Barycentric, kw...)
    f = x -> 0; @test pass(f, approx(f, max_iter=1), pts, atol=tol)
    f = x -> x; @test pass(f, approx(f, max_iter=2), pts, atol=tol)
    f = x -> 1im*x; @test pass(f, approx(f, max_iter=4), pts, atol=tol)
    f = x -> x+x^2; @test pass(f, approx(f, max_iter=3), pts, atol=tol)
    f = x -> x+x^3; @test pass(f, approx(f, max_iter=4), pts, atol=tol)
    f = x -> 1/(1.1+x); @test pass(f, approx(f, max_iter=4), pts, atol=tol)
    f = x -> 1/(1+1im*x); @test pass(f, approx(f, max_iter=4), pts, atol=tol)
    f = x -> 1/(3+x+x^2); @test pass(f, approx(f, max_iter=3), pts, atol=tol)
    f = x -> 1/(1.01+x^3); @test pass(f, approx(f, max_iter=4), pts, atol=tol)
    f = x -> cis(x); r = approx(f);
    @test sum(@. abs(r(pts))-1)/length(pts) < 2e-13
    f = x -> exp(exp(x))/(x - 0.2im); r = approx(f);
    @test minimum(abs.(poles(r) .- .2im)) < 1e-12
end

@testset "Interval [$a, $b]" for (a, b) in ((-2, 3), (0, 4), (-2e-4, 0), (-3e3, 5e6))
    pts = range(a, b, 1000)
    approx(f; kw...) = approximate(f, pts; method=Barycentric, kw...)
    tol(f) = 1000 * eps() * max(b - a, abs(f(a)), abs(f(b)))
    f = x -> 1 / sin(a + (b-a)*1.05im - x); @test pass(f, approx(f), pts, atol=tol(f))
    f = x -> exp(-10/(a + 1.1*(b-a) - x)); @test pass(f, approx(f), pts, atol=tol(f))
end
