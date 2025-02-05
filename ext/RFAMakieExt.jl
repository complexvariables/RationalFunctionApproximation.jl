module RFAMakieExt

using RationalFunctionApproximation, Makie, ComplexRegions, Logging, Printf
RFA = RationalFunctionApproximation

"""
    convergenceplot(r)

Plot the convergence history of a rational approximation.

Markers show the maximum error on (the boundary of) the domain as a function of the numerator/denominator degree. A red marker indicates that the approximation has disallowed poles in its domain. A gold halo highlights the best approximation.
"""
function RFA.convergenceplot(r::RFA.Approximation)
    ismissing(r.stats) && error("No convergence data")
    fig = Figure( )
    ax = Axis(
        fig[1,1], xlabel="degree", ylabel="max error", yscale=log10
        )
    stats = r.stats
    deg = [degree(rewind(r, n)) for n in eachindex(stats.error)]
    color = [any(b) ? :red : :darkblue for b in stats.isbad]
    scatter!(ax, deg, stats.error; color, colorrange=(0, 1), markersize=10)
    scatter!(ax, deg[stats.bestidx], stats.error[stats.bestidx],
        color=RGBAf(1,1,1,0), strokecolor=:gold, strokewidth=3, markersize=16 )
    return fig
end

"""
    errorplot(r; use_abs=false)

Plot the pointwise error of an `Approximation` on (the boundary of) its domain. If the error is not real, then the real and imaginary parts are plotted separately, unless `use_abs=true`.
"""
function RFA.errorplot(r::RFA.Approximation; use_abs=false)
    fig = Figure( )
    ax = Axis(
        fig[1,1], xlabel="boundary parameter", ylabel="error"
        )
    p = r.domain isa RFA.ComplexSCRegion ? r.domain.boundary : r.domain
    t, τ = RFA.refine(p, r.prenodes, max(30, ceil(Int,600/max(degree(r)...))))
    idx = sortperm(t)
    t = t[idx]; τ = τ[idx];
    err = @. r.original(τ) - r.fun.(τ)
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

function RFA.animate(r::RFA.Approximation, filename=tempname()*".mp4")
    function dilate!(box, factor)
        center = (box[2] + box[1]) / 2
        radius = (box[2] - box[1]) / 2
        box[1] = center - factor * radius
        box[2] = center + factor * radius
        return box
    end
    function get_error(n)
        x, y = with_logger(SimpleLogger(stderr, Logging.Warn)) do
            check(rewind(r, n))
        end
        idx = sortperm(x)
        return x[idx], y[idx]
    end
    function max_error_round(err)
        M = maximum(abs, err)
        return round(maximum(abs.(err)), digits=3)
    end
    fig = Figure(size=(850,400) )
    iter = Observable(1)
    ax1 = Axis(fig[1,1], xlabel="Re(z)", ylabel="Im(z)", title=@lift("Degree $($iter-1)"))
    lines!(ax1, discretize(r.domain))
    poles = @lift Point2.(real(r.stats.poles[$iter]), imag(r.stats.poles[$iter]))
    xbox, ybox = ComplexRegions.enclosing_box(discretize(r.domain))
    if xbox[1] == xbox[2]
        xbox = ybox
    end
    if ybox[1] == ybox[2]
        ybox = xbox
    end
    dilate!(xbox, 1.2)
    dilate!(ybox, 1.2)
    nodecolor = @lift([bad ? :red : :black for bad in r.stats.isbad[$iter]])
    scatter!(ax1, poles, color=nodecolor, markersize=7)
    limits!(ax1, xbox..., ybox...)
    err_data = @lift(Point2.(get_error($iter)...))
    max_err = @lift(maximum(abs(e) for e in getindex.($err_data, 2)))
    ax2 = Axis(fig[1,2],
        xlabel="parameterization", ylabel="error",
        title = @lift(@sprintf("max error = %.2e", $max_err))
    )
    max_err_hi = @lift(10 ^ min(0, 2 * ceil(log10($max_err)/2)))
    max_err_lo = @lift(-$max_err_hi)
    lines!(ax2, err_data)
    record(fig, filename, eachindex(r.stats.poles); framerate=2)  do n
        iter[] = n
        ylims!(ax2, max_err_lo[], max_err_hi[])
        # autolimits!(ax2)
    end
    return filename
end

end # module
