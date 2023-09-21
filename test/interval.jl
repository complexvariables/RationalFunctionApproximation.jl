pts = 10 .^ range(-15,0,500); pts = [-reverse(pts); 0; pts]
approx(f; kw...) = approximate(f, unit_interval; kw...)

@testset "Basic functions" begin
    f = x -> abs(x + 0.5 + 0.01im); @test pass(f, approx(f), pts, atol=5e-13)
    f = x -> sin(1/(1.05-x)); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> exp(-1/(x^2)); @test pass(f, approx(f), pts, rtol=4e-13)
    f = x -> exp(-100x^2); @test pass(f, approx(f), pts, rtol=2e-13)
    f = x -> exp(-10/(1.2-x)); @test pass(f, approx(f), pts, rtol=1e-12)
    f = x -> 1/(1+exp(100*(x+.5))); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> sin(100x) * exp(-10x^2); @test pass(f, approx(f), pts, atol=1e-11)
    f = x -> abs(x);  @test pass(f, approx(f), pts, atol=1e-8)
    f = x -> abs(x - 0.95);  @test pass(f, approx(f), pts, atol=1e-6)
end

@testset "Low-accuracy" begin
    f = x -> exp(3x);
    @test !pass(f, approximate(f, unit_interval, tol=1e-4), pts, atol=1e-8)
end

@testset "Poles, zeros, residues" begin
    f = z -> (z+1) * (z+2) / ((z+3) * (z+4))
    r = approx(f)
    pol = poles(r)
    zer = roots(r)
    @test isapprox(sum(pol+zer), -10, atol=1e-12)

    f = x -> 2/(3 + x) + 5/(x - 2im);  r = approx(f)
    @test isapprox(prod(values(residues(r))), 10, atol=1e-8)

    f = x -> sinpi(10x);  r = approx(f);
    @test isapprox(sort(abs.(roots(r)))[19], 0.9, atol=1e-12)

    f = z -> (z - (3 + 3im))/(z + 2);  r = approx(f)
    pol,zer = poles(r), roots(r)
    @test isapprox(pol[1]*zer[1], -6-6im, atol=1e-12)
end

@testset "Vertical scaling" begin
    f = x -> 1e100*sin(x); @test pass(f, approx(f), pts, rtol=2e-13)
    f = x -> 1e-100*cos(x); @test pass(f, approx(f), pts, rtol=2e-13)
    # f = x -> 1e100*sin(x); @test pass(f, approx(f, =5, lawson=20), pts, rtol=1e-6)
    # f = x -> 1e-100*cos(x); @test pass(f, approx(f, =5, lawson=20), pts, rtol=1e-6)
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
    f = x -> tanh(100x); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> tanh(100*(x-.2)); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> exp(x); @test pass(f, approx(f, tol=1e-13), pts, atol=2e-13)
    f = x -> cis(x); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> cis(x); r = approx(f);
    @test sum(@. abs(r(pts))-1)/length(pts) < 2e-13
    f = x -> exp(exp(x))/(x - 0.2im); r = approx(f);
    @test minimum(abs.(poles(r) .- .2im)) < 1e-12
end

@testset "Other intervals" begin
    zz(a,b) = (pts .+ 1)*(b-a)/2 .+ a
    for (a,b) in ((-2,3), (0,1), (-0.01,0) , (-1e4, 2e6))
        for f in ( x -> abs(x + 0.5 + 0.01im), x -> sin(1/(1.05im-x)), x -> exp(-10/(1.2*(b+1)-x)) )
            @test pass(f, approximate(f, Segment(a,b)), zz(a,b), rtol=1e-9)
        end
    end
end
