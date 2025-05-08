module RFAMakieExt

using RationalFunctionApproximation, Makie, ComplexRegions, Logging, Printf
RFA = RationalFunctionApproximation

function RFA.convergenceplot(r::RFA.Approximation)
    deg, err, _, allowed, best = get_history(r)
    fig = Figure( )
    ax = Axis(fig[1,1], xlabel="degree", ylabel="relative max error", yscale=log10)
    color = [all(b) ? :darkblue : :red for b in allowed]
    scatter!(ax, deg, err; color, colorrange=(0, 1), markersize=10)
    scatter!(ax, deg[best], err[best],
        color=RGBAf(1,1,1,0), strokecolor=:gold, strokewidth=3, markersize=16 )
    return fig
end

function RFA.errorplot(r::RFA.Approximation; use_abs=false)
    fig = Figure( )
    ax = Axis(fig[1, 1], xlabel="boundary parameter", ylabel="error")
    t, τ, err = check(r, quiet=true, prenodes=true)
    if use_abs
        lines!(ax, t, abs.(err))
        ax.ylabel = "| error |"
    elseif isreal(err)
        lines!(ax, t, err)
    else
        lines!(ax, t, real.(err))
        lines!(ax, t, imag.(err))
        ax.ylabel = "Re, Im error"
    end
    return fig
end

function axisbox(z)
    xbox, ybox = ComplexRegions.enclosing_box(z, 1.42)
    center = [(box[2] + box[1]) / 2 for box in (xbox, ybox)]
    radius = maximum((box[2] - box[1]) / 2 for box in (xbox, ybox))
    xbox = center[1] .+ [-radius, radius]
    ybox = center[2] .+ [-radius, radius]
    return xbox, ybox
end

function RFA.poleplot(r::RFA.Approximation, idx::Integer=0)
    fig = Figure( )
    ax = Axis(fig[1,1], xlabel="Re(z)", ylabel="Im(z)", aspect=DataAspect())
    z, _ = check(r, quiet=true)
    z = [z; z[1]]
    lines!(ax, real(z), imag(z))
    zp = iszero(idx) ? RFA.poles(r) : RFA.poles(rewind(r, idx))
    color = [r.allowed(z) ? :black : :red for z in zp]
    scatter!(ax, Point2.(real(zp), imag(zp)); color, markersize=7)
    xbox, ybox = axisbox(z)
    limits!(ax, xbox..., ybox...)
    return fig
end

function RFA.animate(r::RFA.Approximation, filename=tempname()*".mp4")
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

    record(fig, filename, eachindex(r.history.len); framerate=2) do n
        iter[] = n
        ylims!(ax2, max_err_lo[], max_err_hi[])
        # autolimits!(ax2)
    end
    return filename
end

end # module
