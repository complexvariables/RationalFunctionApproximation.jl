pts = cispi.(2 * (0:500) / 500)
UD = unit_disk
UC = unit_circle

@testset "Unit disk" begin
    f = z -> sin(10z) * exp(-z^2); @test pass(f, approximate(f, UD), pts, rtol=2e-11)
    f = z -> sin(1/(1.1 - z)); @test pass(f, approximate(f, UD), pts, rtol=2e-13)
    f = sec; @test pass(f, approximate(f, UD, degree=6), pts, rtol=1e-6)
    f = z -> cos(sin(z)) + exp(z)/(z-1.1); @test pass(f, approximate(f,UD), pts, rtol=2e-13)
    f = x -> cis(x);  @test pass(f, approximate(f, UD), pts, atol=2e-13)
end

@testset "Unit circle" begin
    f = z -> abs(z-1im);  @test pass(f, approximate(f, UC), pts, rtol=2e-12)
    f = z -> tan(π*z);  @test pass(f, approximate(f, UC), pts, rtol=2e-13)
    f = z -> tanh(100z); @test pass(f, approximate(f, UC), pts, rtol=2e-13)
end

@testset "Translate and scale" begin
    f = z -> sin(10z) * exp(-z^2)
    for (a,c) in ( (2.5, 0), (1, -1im), (0.4, -2))
        F = approximate(f, a*UC + c)
        @test pass(f, F, a*pts .+ c)
    end
end

@testset "Poles, zeros, residues" begin
    f = z -> tan(π*z);  F = approximate(f, UC)
    pol = poles(F); @test sort(abs.(pol))[1:5] ≈ 0.5*[1;1;3;3;5] atol=1e-3

    f = z -> exp(exp(z)) / (z - 0.2im); pol = poles(approximate(f, UC));
    @test minimum(@. abs(pol - .2im)) < 1e-10

    f = z -> (z+1) * (z+2) / ((z+3) * (z+4));  F = approximate(f, UC)
    pol = poles(F);  zer = roots(F);
    @test isapprox(sum(pol+zer), -10, atol=1e-12)

    f = z -> 2/(3+z) + 5/(z-2im);  F = approximate(f, UD)
    @test isapprox( prod(values(residues(F))), 10, atol=1e-8 )

    f = z -> (z-(3+3im))/(z+2);  F = approximate(f, UD)
    pol, zer = poles(F), roots(F);  @test isapprox(pol[1]*zer[1], -6-6im, atol=1e-12)
end

@testset "Tolerance" begin
    f = z -> exp(3*z); @test !pass(f, approximate(f, UD, tol=1e-4), pts, atol=1e-8)
    f = z -> exp(3*z); @test pass(f, approximate(f, UD, tol=1e-10), pts, atol=1e-8)
end

@testset "Vertical scaling" begin
    f = z -> 1e100sin(z); @test pass(f, approximate(f, UD), pts, rtol=2e-13)
    @test pass(f, approximate(f, UD, degree=5), pts, rtol=1e-6)
    # f = z -> 1e-100cos(z); @test pass(f, approximate(f, UD), pts, rtol=2e-13)
    # @test_broken pass(f, approximate(f, UD, degree=5), pts, rtol=1e-6)
end

# @testset "Lawson" begin
#     f = z -> exp(z); r = @test pass(f, approximate(f, UD, degree=3, lawson=20), pts, atol=1e-3)
#     f = z -> cis(3z); @test pass(f, approximate(f, UD, degree=6, lawson=20), pts, atol=1e-3)
# end

@testset "Polynomials and reciprocals" begin
    for R in (UC, UD)
        for f in (z->0, z->z, z->1im*z, z->z+z^2, z->z+z^3)
            @test pass(f, approximate(f, R), pts, atol=2e-13)
            @test pass(f, approximate(f, R, degree=3), pts, rtol=2e-13)
        end
    end
    for f in (z->1/(1.1+z), z->1/(2+1im*z), z->1/(3+z+z^2), z->1/(1.01+z^3))
        @test pass(f, approximate(f, UC), pts, rtol=2e-13)
    end
end

@testset "Specified degree" begin
    f = x -> x+x^2; @test pass(f, approximate(f, degree=2, UD), pts, atol=2e-13)
    f = x -> x+x^3; @test pass(f, approximate(f, degree=3, UD), pts, atol=2e-13)
    f = x -> x+x^3; @test !pass(f, approximate(f, degree=2, UD), pts, atol=2e-13)
    f = x -> 1/(3im + x + x^2); @test pass(f, approximate(f, degree=2, UC), pts, rtol=2e-13)
    f = x -> 1/(3im + x + x^2); @test !pass(f, approximate(f, degree=1, UC), pts, rtol=2e-13)
    f = x -> 1/(1.01 + x^3); @test pass(f, approximate(f, degree=3, UD), pts, rtol=2e-13)
    f = x -> 1/(1.01 + x^3); @test !pass(f, approximate(f, degree=2, UD), pts, rtol=2e-13)
end
