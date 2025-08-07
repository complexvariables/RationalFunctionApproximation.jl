@testset "Custom path $ipath" for (ipath, c) in enumerate( [
    ClosedCurve(t -> cispi(2t) + 0.2cospi(6t) - 0.1sinpi(4t)),
    Shapes.star,
    cispi(0.31) * Shapes.square + 0.1im,
])
    @test approximate(tan, c) ≈ tan
    @test approximate(cos, interior(c)) ≈ cos
    f = z -> exp(1/z)
    @test approximate(f, exterior(c)) ≈ f
end
