test_points = Dict()
for T in (Float64, Double64)
    z = 1im * T(10) .^ range(T(-15), T(0), 500);
    test_points[T] = [-reverse(z); 0; z]
end
domain = Dict()
for T in (Float64, Double64)
    domain[T] = Segment{T}(-1im, 1im)
end

@testset "Basic functions for $method" for method in (Barycentric, Thiele)
    T = Float64
    tol = 2000*eps(T)
    pts = test_points[T]
    approx(f; kw...) = approximate(f, domain[T]; method, kw...)
    f = x -> abs(x + 1//2 + 1im//100); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> sinh(1 / (21//20 - x)); @test pass(f, approx(f), pts; rtol=tol)
    if method != Thiele
        f = x -> x + exp(-1 / x^2); @test pass(f, approx(f; stagnation=30), pts; rtol=tol) skip=true
        f = x -> x + tan(100x); @test pass(f, approx(f), pts, rtol=tol)
    end
    f = x -> x + exp(100x^2); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> exp(-10 / (6//5 - x)); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> x + sinh(80x) * exp(10x^2); @test pass(f, approx(f; stagnation=30), pts; rtol=tol)
    f = x -> 10x + tan(100*(x - 1//5)); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> exp(x); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> cis(x); @test pass(f, approx(f), pts, rtol=tol)
end

@testset "Double64 for $method" for method in (Barycentric, Thiele)
    T = Double64
    pts = test_points[T]
    tol = 2000*eps(T)
    approx(f; kw...) = approximate(f, domain[T]; method, kw...)
    f = x -> abs(x + 1//2 + 1im//100); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> sinh(1 / (21//20 - x)); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> exp(-100*(x + 1//4)^2); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> exp(-10 / (6//5 - x)); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> 1/(1 + exp(100*(x + 1//2))); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> x + sinh(50x) * exp(10x^2); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> x + tan(100*(x - 1//5)); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> exp(x); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> cis(x); @test pass(f, approx(f), pts, rtol=tol)
end

@testset "Low accuracy" begin
    T = Float64
    method = Barycentric
    pts = test_points[T]
    approx(f; kw...) = approximate(f, domain[T]; method, kw...)
    f = x -> exp(3x);
    r = approx(f, tol=1e-5)
    @test !pass(f, r, pts, atol=1e-10)
    @test pass(f, r, pts, atol=5e-5)
    f = x -> abs(x);  @test pass(f, approx(f, stagnation=30), pts, atol=1e-10)
    f = x -> abs(x - 0.95);  @test pass(f, approx(f, stagnation=30), pts, atol=1e-9)
end

@testset "Poles, zeros, residues" for T in (Float64,)
    approx(f; kw...) = approximate(f, domain[T]; kw...)
    f = z -> (z+1) * (z+2) / ((z+3) * (z+4))
    r = approx(f)
    pol = poles(r)
    zer = roots(r)
    @test isapprox(sum(pol+zer), -10, rtol=5000*eps(T))

    f = x -> 2/(3 + x) + 5/(x - 2im);  r = approx(f)
    @test isapprox(prod(residues(r)[2]), 10, rtol=sqrt(eps(T)))

    f = x -> sinh(10pi*x);  r = approx(f);
    zer = roots(r)
    miss = minimum(abs, zer .- 9im//10)
    @test miss < 1e6*eps(T)

    f = z -> (z - (3 + 3im))/(z + 2);  r = approx(f)
    pol, zer = poles(r), roots(r)
    @test isapprox(pol[1]*zer[1], -6-6im, rtol=5000*eps(T))
end

@testset "Vertical scaling in $T" for T in (Float64, Double64)
    approx(f; kw...) = approximate(f, domain[T]; method=Barycentric, kw...)
    pts = test_points[T]
    f = x -> T(10)^50*sinh(x); @test pass(f, approx(f), pts, rtol=2000*eps(T))
    f = x -> T(10)^(-50)*cosh(x); @test pass(f, approx(f), pts, rtol=2000*eps(T))
end

# @testset "Lawson" begin
#     f = x -> exp(x); @test pass(f, aaax(f, max_iter=3, lawson=20), xx, atol=1e-3)
#     f = x -> cis(3x); @test pass(f, aaax(f, max_iter=3, lawson=20), xx, atol=1e-3)
#     f = z -> exp(z); r = @test pass(f, aaaz(f, max_iter=3, lawson=20), zz, atol=1e-3)
#     f = z -> cis(3z); @test pass(f, aaaz(f, max_iter=6, lawson=20), zz, atol=1e-3)
# end

@testset "Polynomials and reciprocals" begin
    T = Float64
    pts = test_points[T]
    tol = 2000*eps(T)
    approx(f; kw...) = approximate(f, domain[T]; method=Barycentric, kw...)
    f = x -> 0; @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> x; @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> 1im*x; @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> x + x^3; @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> x + x^2; @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> 1/(1.1 + x); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> 1/(1 + x); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> 1/(3 + x + x^2); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> 1/(1.01 + x^3); @test pass(f, approx(f), pts, atol=2e-13)
end
