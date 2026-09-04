# The pre-rename API selected the interpolant type with a `method` keyword argument.
# That spelling still works, but must announce itself: `Base.depwarn` is silent under
# Julia's default `--depwarn=no`, so the shim passes `force=true`.
@testset "Deprecated `method` keyword" begin
    f = exp
    z = collect(range(-1, 1, 200))
    warned = (:warn, r"deprecated")

    @testset "warns on every entry point" begin
        @test_logs warned match_mode=:any approximate(f, Segment(-1, 1); method=Barycentric)
        @test_logs warned match_mode=:any approximate(f, -1..1; method=Barycentric)
        @test_logs warned match_mode=:any approximate(f, interior(Circle(0, 1)); method=Barycentric)
        @test_logs warned match_mode=:any approximate(f, z; method=Barycentric)
        @test_logs warned match_mode=:any approximate(f.(z), z; method=Barycentric)
    end

    @testset "selects the requested type" begin
        @test approximate(f, Segment(-1, 1); method=Barycentric).fun isa Barycentric
        @test approximate(f, Segment(-1, 1); method=Thiele).fun isa Thiele
        # an instance is accepted as well as a type
        @test approximate(f, Segment(-1, 1); method=Barycentric()).fun isa Barycentric
    end

    @testset "other keywords survive the shim" begin
        r = approximate(x -> exp(3x), Segment(-1, 1); method=Barycentric, tol=1e-5)
        @test r.fun isa Barycentric
        @test !pass(x -> exp(3x), r, range(-1, 1, 500), atol=1e-7)
        @test pass(x -> exp(3x), r, range(-1, 1, 500), atol=5e-5)
    end

    @testset "no warning without the keyword" begin
        @test_logs approximate(f, Segment(-1, 1))
        @test_logs approximate(f, Segment(-1, 1), Barycentric())
    end
end
