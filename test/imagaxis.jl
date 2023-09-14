@testset "varied and random stuff" begin
    _,_,_,_,err,_ = aaax( x->abs(x) ) 
    @test abs(err) < 1e-3
    _,pol,_ = aaaz( z->cos(sin(z)) + exp(z)/(z-11//10) )
    @test isapprox(minimum(abs.(pol)), 1.1, atol=1e-8)
    f = x -> abs(x + 1//2 + 1im//100)
    @test pass(f, aaax(f), xx, atol=2e-13)
end

@testset "Poles and zeros" begin
    f = z -> (z+1) * (z+2) / ((z+3) * (z+4))
    for aaa in active
        r,pol,res,zer,err,_ = aaa(f)
        @test isapprox(sum(pol+zer), -10, atol=1e-12)
    end

    f = z -> 1 / (sqrt(z+1) * sqrt(z+2)); @test pass(f, aaai(f), ii, atol=2e-13)
    f = z -> 11 / (sqrt(z + 10im + 1//10) * sqrt(z - 10im + 1//10))
    @test pass(f, aaai(f), ii, atol=2e-12)
end

@testset "Low-accuracy computations" begin
    f = x -> exp(3*x); @test !pass(f, aaax(f, degree=150, lawson=0, tol=1e-4), xx, atol=1e-8)
    f = z -> exp(3*z); @test !pass(f, aaaz(f, degree=150, lawson=0, tol=1e-4), zz, atol=1e-8)
    f = z -> 1 / sqrt(z+3); @test !pass(f, aaai(f, degree=150, lawson=0, tol=1e-4), ii, atol=1e-8)
end

@testset "Miscellaneous functions" begin
    f = z -> sin(10z) * exp(-z^2); @test pass(f, aaaz(f), zz, rtol=2e-11) 
    f = z -> sin(1/(1.1 - z)); @test pass(f, aaaz(f), zz, rtol=2e-13)
    f = sec; @test pass(f, aaaz(f,degree=6), zz, rtol=1e-6) 
    f = z -> z^4; @test pass(f, aaaz(f,degree=3), zz, rtol=20) 
    f = z -> tan(π*z); r,pol,_ = aaaz(f,degree=150,lawson=0,tol=1e-13,mero=true)
    @test sort(abs.(pol))[1:5] ≈ 0.5*[1;1;3;3;5] atol=1e-3

    f = x -> sin(1/(1.05-x)); @test pass(f, aaax(f), xx, atol=2e-13) 
    f = x -> exp(-1/(x^2)); @test pass(f, aaax(f), xx, rtol=2e-13)
    f = x -> exp(-100*x^2); @test pass(f, aaax(f), xx, rtol=2e-13)
    f = x -> exp(-10/(1.2-x)); @test pass(f, aaax(f), xx, rtol=1e-12)
    f = x -> exp(-10/(1.2-x)); @test pass(f, aaax(f,degree=8,lawson=20), xx, rtol=1e-8)
    # f = x -> gamma(x+2); [r,pol] = aaax(f); pass(f,abs(min(pol+2))<1e-5)
    f = x -> 1/(1+exp(100*(x+.5))); @test pass(f, aaax(f), xx, atol=2e-13)
end

@testset "Zeros and residues" begin
    f = z -> 2/(3+z) + 5/(z-2im); r,pol,res,_ = aaaz(f)
    @test isapprox( prod(res), 10, atol=1e-8 )
    f = x -> 2/(3+x) + 5/(x-2im); r,pol,res,_ = aaax(f); @test isapprox(prod(res), 10, atol=1e-8)
    f = x -> sinpi(10x); r,pol,res,zer,_ = aaax(f); @test isapprox(sort(abs.(zer))[19], 0.9, atol=1e-12)

    f = z -> (z-(3+3im))/(z+2)
    for aaa in (aaax, aaaz, aaai)
        r,pol,res,zer,_ = aaa(f, degree=150, lawson=0, tol=1e-13) 
        @test isapprox(pol[1]*zer[1], -6-6im, atol=1e-12)
    end
end

@testset "Imaginary axis" begin 
    for scl in (1.0, 1e100, 1e-100)
        f = z -> scl/((z+(1+2im))*(z+.1)); r,pol,_ = aaai(f)
        @test isapprox(prod(pol),0.1+0.2im, atol=1e-11)
        f = z -> scl/((z-(1+2im))*(z+.1)); r,pol,_ = aaai(f, mero=true)
        @test isapprox(prod(pol),-0.1-0.2im, atol=1e-11)
    end
end

@testset "Vertical scaling" begin
    f = x -> 1e100*sin(x); @test pass(f, aaax(f), xx, rtol=2e-13)
    f = x -> 1e-100*cos(x); @test pass(f, aaax(f), xx, rtol=2e-13)
    f = z -> 1e100*sin(z); @test pass(f, aaaz(f), zz, rtol=2e-13)
    f = z -> 1e-100*cos(z); @test pass(f, aaaz(f), zz, rtol=2e-13) 
    f = x -> 1e100*sin(x); @test pass(f, aaax(f, degree=5, lawson=20), xx, rtol=1e-6) 
    f = x -> 1e-100*cos(x); @test pass(f, aaax(f, degree=5, lawson=20), xx, rtol=1e-6) 
    f = z -> 1e100*sin(z); @test pass(f, aaaz(f,degree=5,lawson=20), zz, rtol=1e-6)
    f = z -> 1e-100*cos(z); @test pass(f, aaaz(f,degree=5,lawson=20), zz, rtol=1e-6)
end

@testset "Lawson" begin
    f = x -> exp(x); @test pass(f, aaax(f, degree=3, lawson=20), xx, atol=1e-3)
    f = x -> cis(3x); @test pass(f, aaax(f, degree=3, lawson=20), xx, atol=1e-3)
    f = z -> exp(z); r = @test pass(f, aaaz(f, degree=3, lawson=20), zz, atol=1e-3)
    f = z -> cis(3z); @test pass(f, aaaz(f, degree=6, lawson=20), zz, atol=1e-3)
end

@testset "Polynomials and reciprocals" begin
    args = Dict(:degree=>150, :lawson=>0, :tol=>1e-13)
    f = x -> 0; @test pass(f, aaax(f; args...), xx, atol=2e-13)
    f = x -> x; @test pass(f, aaax(f; args...), xx, atol=2e-13)
    f = x -> 1im*x; @test pass(f, aaax(f; args...), xx, atol=2e-13)
    f = x -> x+x^2; @test pass(f, aaax(f; args...), xx, atol=2e-13)
    f = x -> x+x^3; @test pass(f, aaax(f; args...), xx, atol=2e-13) 
    f = x -> 1/(1.1+x); @test pass(f, aaax(f; args...), xx, atol=2e-13) 
    f = x -> 1/(1+1im*x); @test pass(f, aaax(f; args...), xx, atol=2e-13)
    f = x -> 1/(3+x+x^2); @test pass(f, aaax(f; args...), xx, atol=2e-13)
    f = x -> 1/(1.01+x^3); @test pass(f, aaax(f; args...), xx, atol=2e-13)
    
    args = Dict(:degree=>150, :lawson=>0, :tol=>1e-13, :mero=>true)
    f = z -> 0; @test pass(f, aaaz(f; args...), zz, atol=2e-13)
    f = z -> z; @test pass(f, aaaz(f; args...), zz, atol=2e-13)
    f = z -> 1im*z; @test pass(f, aaaz(f; args...), zz, atol=2e-13)
    f = z -> z+z^2; @test pass(f, aaaz(f; args...), zz, atol=2e-13)
    f = z -> z+z^3; @test pass(f, aaaz(f; args...), zz, atol=2e-13)
    f = z -> 1/(1.1+z); @test pass(f, aaaz(f; args...), zz, atol=2e-13)
    f = z -> 1/(2+1im*z); @test pass(f, aaaz(f; args...), zz, atol=2e-13)
    f = z -> 1/(3+z+z^2); @test pass(f, aaaz(f; args...), zz, atol=2e-13)
    f = z -> 1/(1.01+z^3); @test pass(f, aaaz(f; args...), zz, atol=1e-11)
end

@testset "Specified degree" begin
    f = x -> 0; @test pass(f, aaax(f, degree=0, lawson=0), xx, atol=2e-13) 
    f = x -> x; @test pass(f, aaax(f, degree=1, lawson=0), xx, atol=2e-13) 
    f = x -> 1im*x; @test pass(f, aaax(f, degree=3, lawson=0), xx, atol=2e-13) 
    f = x -> x+x^2; @test pass(f, aaax(f, degree=2, lawson=0), xx, atol=2e-13) 
    f = x -> x+x^3; @test pass(f, aaax(f, degree=3, lawson=0), xx, atol=2e-13) 
    f = x -> 1/(1.1+x); @test pass(f, aaax(f, degree=3, lawson=0), xx, atol=2e-13) 
    f = x -> 1/(1+1im*x); @test pass(f, aaax(f, degree=3, lawson=0), xx, atol=2e-13) 
    f = x -> 1/(3+x+x^2); @test pass(f, aaax(f, degree=2, lawson=0), xx, atol=2e-13) 
    f = x -> 1/(1.01+x^3); @test pass(f, aaax(f, degree=3, lawson=0), xx, atol=2e-13) 
    f = x -> tanh(100x); @test pass(f, aaax(f), xx, atol=2e-13) 
    f = x -> tanh(100*(x-.2)); @test pass(f, aaax(f), xx, atol=2e-13) 
    f = x -> exp(x); @test pass(f, aaax(f, tol=1e-13), xx, atol=2e-13) 
    f = x -> cis(x); @test pass(f, aaax(f), xx, atol=2e-13) 
    f = x -> cis(x); r,_ = aaax(f); 
    @test mean(@. abs(r(xx))-1) < 2e-13
    f = x -> exp(exp(x))/(x-0.2im); r,pol,_ = aaax(f); 
    @test minimum(@. abs(pol-.2im)) < 1e-10
end

@testset "Slower" begin
    f = x -> sin(100x) * exp(-10x^2); @test pass(f, aaax(f), xx, atol=1e-11)
    f = x -> abs(x);  @test pass(f, aaax(f), xx, atol=1e-8)
    f = x -> abs(x - 0.95);  @test pass(f, aaax(f), xx, atol=1e-6)
end

# this set requres DoubleFloats to be installed
# functions must not introduce any standard floats
# (note that // creates exact rational numbers)
@testset "Extended precision" begin 
    using DoubleFloats
    T = Double64
    xx = T(10) .^ range(T(-15),T(0),500); xx = [-reverse(xx); 0; xx]
    zz = cispi.(2 * T.((0:500) / 500))
    for f in [x -> tanh(20x), x -> 1im + abs(x + 1//2 + 1im//10) ]
        @test pass(f, aaax(f, numtype=Double64), xx, atol=8e-29 )
        @test pass(f, aaaz(f, mero=true, numtype=Double64), zz, atol=8e-29 )
    end
end