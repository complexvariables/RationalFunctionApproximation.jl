@testset "Derivatives for $method" for (method, domain) in Iterators.product((Barycentric, Thiele), (unit_interval, unit_disk, Shapes.square))
    for (f, df) in (
        (exp, exp),
        (x -> exp(-x), x -> -exp(-x)),
        (x -> cis(x), x -> 1im * cis(x)),
        (x -> x, x -> 1),
        (x -> 1im * x^2, x -> 2im * x),
        (x -> 1 / (1.1 - x), x -> 1 / (1.1 - x)^2),
        (x -> log(1.1 - x), x -> -1 / (1.1 - x)),
        (sin, cos),
    )
        r = approximate(f, domain; method)
        @test(isapprox(derivative(r; allowed=true), df))
    end
end
