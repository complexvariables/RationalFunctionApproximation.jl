module RFAMakieExt

using RationalFunctionApproximation, Makie, ComplexRegions, Printf
RFA = RationalFunctionApproximation

function RFA.convergenceplot!(target, r::RFA.AbstractApproximation; show_best=true, kwargs...)
    deg, err, _, allowed, best = get_history(r)
    attr = Makie.Attributes(; markersize=8)
    if r.allowed != true
        attr = merge(attr, Makie.Attributes(
            markercolor= [all(b) ? :darkblue : :red for b in allowed],
            color=RGBAf(0, 0, 0, 0.5)
            ))
    end
    attr = merge(attr, Makie.Attributes(; kwargs...))
    out1 = scatterlines!(target, deg, err; attr...)
    out2 = if show_best
        scatterlines!(target, deg[best], err[best],
                color=RGBAf(1,1,1,0), strokecolor=:gold, strokewidth=3, markersize=13 )
    else
        nothing
    end
    return out1, out2
end

function RFA.convergenceplot(pos::GridPosition, r::RFA.AbstractApproximation; kwargs...)
    ax = Axis(pos, xlabel="degree", ylabel="relative max error", yscale=log10)
    return Makie.AxisPlot(ax, RFA.convergenceplot!(ax, r; kwargs...)[1])
end

function RFA.convergenceplot(r::RFA.AbstractApproximation; kwargs...)
    fig = Figure( )
    return Makie.FigureAxisPlot(fig, RFA.convergenceplot(fig[1, 1], r; kwargs...)...)
end

function RFA.errorplot!(target, r::RFA.ContinuumApproximation; use_abs=false, kwargs...)
    # try to get enough points for a smooth result
    N = ceil(Int, 1000 / length(nodes(r)))
    t, τ, err = check(r, refinement=N, quiet=true, prenodes=true)
    out = if use_abs
        lines!(target, t, abs.(err); kwargs...)
    elseif isreal(err)
        lines!(target, t, err; kwargs...)
    else
        series!(target, t, hcat(reim(err)...)'; kwargs...)
    end
    return out
end

function RFA.errorplot(pos::GridPosition, r::RFA.ContinuumApproximation; use_abs=false, kwargs...)
    ax = if use_abs
        Axis(pos, xlabel="boundary parameter", ylabel="| error |")
    elseif isreal(r.fun)
        Axis(pos, xlabel="boundary parameter", ylabel="error")
    else
        Axis(pos, xlabel="boundary parameter", ylabel="Re, Im error")
    end
    return Makie.AxisPlot(ax, RFA.errorplot!(ax, r; use_abs, kwargs...))
end

function RFA.errorplot(r::RFA.ContinuumApproximation; kwargs...)
    fig = Figure( )
    return Makie.FigureAxisPlot(fig, RFA.errorplot(fig[1, 1], r; kwargs...)...)
end

function axisbox(z)
    xbox, ybox = ComplexRegions.enclosing_box(z, 1.42)
    center = [(box[2] + box[1]) / 2 for box in (xbox, ybox)]
    radius = maximum((box[2] - box[1]) / 2 for box in (xbox, ybox))
    xbox = center[1] .+ [-radius, radius]
    ybox = center[2] .+ [-radius, radius]
    return xbox, ybox
end

function RFA.poleplot!(target, r::RFA.AbstractApproximation; kwargs...)
    zp = RFA.poles(r)
    attr = Makie.Attributes(; marker=:+, markersize=7, kwargs...)
    if r.allowed != true
        attr = merge(attr, Makie.Attributes(
            color= [r.allowed(z) ? :black : :red for z in zp],
            ))
    end
    return scatter!(target, Point2.(real(zp), imag(zp)); attr...)
end

function RFA.poleplot(pos::GridPosition, r::RFA.AbstractApproximation; kwargs...)
    ax = Axis(pos, xlabel="Re(z)", ylabel="Im(z)", aspect=DataAspect())
    out = RFA.poleplot!(ax, r; kwargs...)
    xbox, ybox = axisbox(test_points(r))
    limits!(ax, xbox..., ybox...)
    return Makie.AxisPlot(ax, out)
end

function RFA.poleplot(r::RFA.AbstractApproximation; kwargs...)
    fig = Figure( )
    return Makie.FigureAxisPlot(fig, RFA.poleplot(fig[1, 1], r; kwargs...)...)
end

function RFA.animate(r::RFA.AbstractApproximation, filename=tempname()*".mp4")
    function get_error(n, x)
        y = r.original.(x) - rewind(r, n).(x)
        return abs.(y)
    end

    t, x, _ = check(r, quiet=true, prenodes=true)
    xbox, ybox = axisbox(x)

    fig = Figure(size=(850,400) )
    iter = Observable(1)

    # left axis: boundary and poles
    ax1 = Axis(fig[1,1],
        xlabel="Re(z)", ylabel="Im(z)",
        title=@lift("Degree $(degree(rewind(r, $iter)))"),
        aspect=DataAspect()
        )
    lines!(ax1, real(x), imag(x))
    limits!(ax1, xbox..., ybox...)
    zp = @lift complex(RFA.poles(rewind(r, $iter)))
    pole = @lift Point2.(real($zp), imag($zp))
    color = @lift [r.allowed(z) ? :black : :red for z in $zp]
    scatter!(ax1, pole; color, markersize=7)

    # right axis: error
    err_data = @lift(Point2.(t, get_error($iter, x)))
    max_err = @lift(maximum(abs(e) for e in getindex.($err_data, 2)))
    ax2 = Axis(fig[1,2],
        xlabel="parameterization", ylabel="error",
        title = @lift(@sprintf("max error = %.2e", $max_err))
    )
    max_err_hi = @lift(10 ^ min(0, 2 * ceil(log10($max_err)/2)))
    max_err_lo = @lift(-$max_err_hi / 20)
    lines!(ax2, err_data)

    record(fig, filename, eachindex(r.history); framerate=2) do n
        iter[] = n
        ylims!(ax2, max_err_lo[], max_err_hi[])
        # autolimits!(ax2)
    end
    return filename
end

end # module
