
@testset "Basic functions in $T" for T in (Float64, Double64)
    pts = T(10) .^ range(T(-15), T(0), 500); pts = [-reverse(pts); 0; pts]
    tol = 2000*eps(T)
    approx(f; kw...) = approximate(f, Segment{T}(-1, 1); kw...)
    f = x -> abs(x + 1//2 + 1im//100); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> sin(1 / (21//20 - x)); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> exp(-1 / x^2); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> exp(-100x^2); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> exp(-10 / (6//5 - x)); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> 1/(1 + exp(100*(x + 1//2))); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> sin(100x) * exp(-10x^2); @test pass(f, approx(f), pts; rtol=tol)
    f = x -> tanh(100*(x - 1//5)); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> tanh(100x); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> exp(x); @test pass(f, approx(f), pts, rtol=tol)
    f = x -> cis(x); @test pass(f, approx(f), pts, rtol=tol)
end

@testset "Low accuracy" begin
    pts = 10.0 .^ range(-15, 0, 500); pts = [-reverse(pts); 0; pts]
    approx(f; kw...) = approximate(f, unit_interval; kw...)
    f = x -> exp(3x);
    r = approximate(f, unit_interval, tol=1e-5)
    @test !pass(f, r, pts, atol=1e-7)
    @test pass(f, r, pts, atol=5e-5)
    f = x -> abs(x);  @test pass(f, approx(f), pts, atol=1e-8)
    f = x -> abs(x - 0.95);  @test pass(f, approx(f), pts, atol=1e-6)
end

@testset "Poles, zeros, residues" for T in (Float64,)
    setprecision(80)
    approx(f; kw...) = approximate(f, Segment{T}(-1, 1); kw...)
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
    pts = T(10) .^ range(T(-15), T(0), 500); pts = [-reverse(pts); 0; pts]
    f = x -> T(10)^50*sin(x); @test pass(f, approx(f), pts, rtol=2000*eps(T))
    f = x -> T(10)^(-50)*cos(x); @test pass(f, approx(f), pts, rtol=2000*eps(T))
end

# @testset "Lawson" begin
#     f = x -> exp(-10/(1.2-x)); @test pass(f, approx(f, =8, lawson=20), pts, rtol=1e-8)
#     f = x -> exp(x); @test pass(f, approx(f, =3, lawson=20), pts, atol=1e-3)
#     f = x -> cis(3x); @test pass(f, approx(f, =3, lawson=20), pts, atol=1e-3)
# end

@testset "Polynomials and reciprocals" begin
    args = Dict(:max_degree=>150, :tol=>1e-13)
    f = x -> 0; @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> x; @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1im*x; @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> x + x^2; @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> x + x^3; @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1/(1.1 + x); @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1/(1 + 1im*x); @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1/(3 + x + x^2); @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1/(1.01 + x^3); @test pass(f, approx(f; args...), pts, atol=2e-13)
end

@testset "Specified degree" begin
    f = x -> 0; @test pass(f, approx(f, max_degree=0), pts, atol=2e-13)
    f = x -> x; @test pass(f, approx(f, max_degree=1), pts, atol=2e-13)
    f = x -> 1im*x; @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
    f = x -> x+x^2; @test pass(f, approx(f, max_degree=2), pts, atol=2e-13)
    f = x -> x+x^3; @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
    f = x -> 1/(1.1+x); @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
    f = x -> 1/(1+1im*x); @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
    f = x -> 1/(3+x+x^2); @test pass(f, approx(f, max_degree=2), pts, atol=2e-13)
    f = x -> 1/(1.01+x^3); @test pass(f, approx(f, max_degree=3), pts, atol=2e-13)
    f = x -> cis(x); r = approx(f);
    @test sum(@. abs(r(pts))-1)/length(pts) < 2e-13
    f = x -> exp(exp(x))/(x - 0.2im); r = approx(f);
    @test minimum(abs.(poles(r) .- .2im)) < 1e-12
end

@testset "Other intervals" begin
    pts = 10 .^ range(-15, 0, 500); pts = [-reverse(pts); 0; pts]
    zz(a,b) = (pts .+ 1)*(b-a)/2 .+ a
    for (a,b) in ((-2,3), (0,1), (-0.01,0) , (-1e4, 2e6))
        for f in ( x -> abs(x + 0.5 + 0.01im), x -> sin(1/(1.05im-x)), x -> exp(-10/(1.2*(b+1)-x)) )
            @test pass(f, approximate(f, Segment(a,b)), zz(a,b), rtol=1e-9)
        end
    end
end
