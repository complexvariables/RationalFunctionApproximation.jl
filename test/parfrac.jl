@testset "Accessor functions" begin
    f, ζ = (x -> 1im*tan(x), [-π/2, π/2, -3π/2, 3π/2]);
    r = approximate(f, Segment(-1, 1), ζ; degree=30)
    @test isapprox(poles(r), [-π/2, π/2, -3π/2, 3π/2])
    @test length(residues(r)[2]) == 4
end


@testset "Real axis" begin
    T = Float64
    tol = 1e-6
    domain[T] = Segment{T}(-1, 1)
    function points(T=Float64)
        z = T(10) .^ range(T(-15), T(0), 500);
        return [-reverse(z); 0; z]
    end
    pts = points(T)
    approx(f, ζ; kw...) = approximate(f, domain[T], ζ; degree=30, kw...)

    f, ζ = (x -> 1im*tan(x), [-π/2, π/2, -3π/2, 3π/2]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> tanh(x), 1im*[-π/2, π/2, -3π/2, 3π/2]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> 1 / sin(21//20 - x), [21//20]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> exp(-10*(x-0.4)^2), []); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> exp(-10 / (6//5 - x)), [6//5]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> exp(x), []); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> exp(x), []); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> 1 / sin(x - 789//788), [789//788]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> x + 2im*x^2, []); @test pass(f, approx(f, ζ; degree=3), pts, atol=2e-13)
    f, ζ = (x -> 1 / (101//100 + x), [-101//100]); @test pass(f, approx(f, ζ; degree=4), pts; rtol=1e-12)
end

@testset "Imaginary axis" begin
    T = Float64
    domain[T] = Segment{T}(-1im, 1im)
    function points(T=Float64)
        z = T(10) .^ range(T(-15), T(0), 500);
        return 1im*[-reverse(z); 0; z]
    end
    tol = 1e-6
    pts = points(T)
    approx(f, ζ; kw...) = approximate(f, domain[T], ζ; degree=30, kw...)

    f, ζ = (x -> 1im*tan(x), [-π/2, π/2, -3π/2, 3π/2]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> tanh(x), 1im*[-π/2, π/2, -3π/2, 3π/2]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> 1 / sin(21//20 - x), [21//20]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> exp(-10*(x-0.4)^2), []); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> exp(-10 / (6//5 - x)), [6//5]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> exp(x), []); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> cis(x), []); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> x + x^2, []); @test pass(f, approx(f, ζ), pts, atol=2e-13)
end

@testset "Double64" begin
    T = Double64
    tol = 1e-18
    domain[T] = Segment{T}(-1, 1)
    function points(T=Float64)
        z = T(10) .^ range(T(-15), T(0), 500);
        return [-reverse(z); 0; z]
    end
    pts = points(T)
    approx(f, ζ; kw...) = approximate(f, domain[T], ζ; degree=80, kw...)

    f, ζ = (x -> tanh(x), 1im*[-π/2, π/2, -3π/2, 3π/2]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> 1 / sin(21//20 - x), [21//20]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> exp(-10 / (6//5 - x)), [6//5]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> cis(x), []); @test pass(f, approx(f, ζ), pts; rtol=tol)
end

@testset "Unit circle" begin
    T = Float64
    tol = 1e-6
    domain[T] = Circle{T}(0, 1)
    pts = [cispi(T(2t)) for t in 0:1//500:1]
    approx(f, ζ; kw...) = approximate(f, domain[T], ζ; degree=30, kw...)

    f, ζ = (x -> 1im*tan(x), [-π/2, π/2, -3π/2, 3π/2]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> tanh(x), 1im*[-π/2, π/2, -3π/2, 3π/2]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> 1 / sin(21//20 - x), [21//20]); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> exp(x), []); @test pass(f, approx(f, ζ), pts; rtol=tol)
    f, ζ = (x -> cis(x), []); @test pass(f, approx(f, ζ), pts; rtol=tol)
end

@testset "Other intervals" begin
    for (a,b) in ((-3,4), (0,1), (-0.01,0) )
        dom = Segment(a,b)
        pts = range(a, b, 800)
        approx(f, ζ; kw...) = approximate(f, dom, ζ; degree=30, kw...)
        f, ζ = (x -> tanh(x), 1im*[-π/2, π/2, -3π/2, 3π/2, -5π/2, 5π/2]);
        @test pass(f, approx(f, ζ), pts; rtol=1e-9)
    end
end

@testset "Arnoldi polynomial" begin
     for path in [Segment(-0.3,1), Circle(1, 2), Shapes.square]
        z = point(path, 0:0.001:1)
        B = RFA.ArnoldiBasis(z, 30)
        @test maximum(abs, (B\cis).(z) - cis.(z)) < 1e-10
     end
end

@testset "Discrete domains" begin
    tol = 1e-6
    approx(f, pts, ζ; kw...) = approximate(f, pts, ζ; degree=30, kw...)
    for path in [Segment(-0.3,1), Circle(1, 2), Shapes.square]
        pts = point(path, 0:0.001:1)
        f, ζ = (x -> 1im*tan(x), [-π/2, π/2, -3π/2, 3π/2]); @test pass(f, approx(f, pts, ζ), pts; rtol=tol)
        f, ζ = (x -> tanh(x), 1im*[-π/2, π/2, -3π/2, 3π/2]); @test pass(f, approx(f, pts, ζ), pts; rtol=tol)
        f, ζ = (x -> exp(x), []); @test pass(f, approx(f, pts, ζ), pts; rtol=tol)
        f, ζ = (x -> cis(x), []); @test pass(f, approx(f, pts, ζ), pts; rtol=tol)
    end
end
