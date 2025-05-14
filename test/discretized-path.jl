@testset "Interval in $T" for T in (Float64, Rational)
    domain = Segment(T(0), T(1))
    p = DiscretizedPath(domain, 5)
    t, z = collect(p, :nodes)
    @test length(t) == length(z) == 5
    @test isapprox(z, (0:4) .// 4)
    @test isapprox(t, (0:4) .// 4)
    @test_throws BoundsError RFA.add_node!(p, [1, 2])

    p = DiscretizedPath(domain, (0:4) .// 4; refinement=4)
    @test size(p.params) == (5, 5)
    @test size(p.points) == (5, 5)
    @test p.next == [2, 3, 4, 0]
    @test p.params[1:4, 3] ≈ @. (0:3) // 4 + (2//20)

    A = reshape(T.(0:19) / 20, 5, 4)
    t, z = collect(p, :all)
    @test isapprox(t[1:end-1], vec(A))
    t, z = collect(p, :test)
    @test isapprox(t, vec(A[2:end, :]))

    domain = Segment(-3, 2)
    p = DiscretizedPath(domain, 5)
    t, z = collect(p, :nodes)
    @test length(t) == length(z) == 5
    @test isapprox(z, @. -3 + 5*(0:4) // 4)
    @test isapprox(t, (0:4) .// 4)
end

@testset "Adding points" begin
    domain = Segment(0.0, 1.0)
    p = DiscretizedPath(domain, (0:4) .// 4; refinement=4, maxpoints=100)
    idx = RFA.add_node!(p, [1, 4])
    @test idx[1:1, :] == collect(CartesianIndices((1:1, 1:5)))
    @test idx[2:2, :] == collect(CartesianIndices((5:5, 1:5)))
    @test p.next == [5, 3, 4, 0, 2]
    @test p.params[1, :] ≈ @. (0:4) * 3 // 100
    @test p.params[5, :] ≈ @. (0:4) // 50 + 3 // 20
    t, z = collect(p)
    @test t ≈ [0, 3, 5, 10, 15, 20] // 20

    idx = RFA.add_node!(p, [5, 1])
    t, z = collect(p)
    @test t ≈ [0, 3, 5, 10, 15, 20] // 20

    idx = RFA.add_node!(p, [5, 5])
    t, z = collect(p)
    @test t ≈ [0, 15, 23, 25, 50, 75, 100] // 100
end
