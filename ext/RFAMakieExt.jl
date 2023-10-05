module RFAMakieExt

using RationalFunctionApproximation, Makie
RFA = RationalFunctionApproximation

"""
    convergenceplot(r)

Plot the convergence history of a `Barycentric` or `Approximation` rational function.

Markers show the maximum error on (the boundary of) the domain as a function of the numerator/denominator degree. A red marker indicates that the approximation has disallowed poles in its domain. A gold halo highlights the best approximation.
"""
RFA.convergenceplot(r::RFA.Approximation) = RFA.convergenceplot(r.fun)
function RFA.convergenceplot(r::Barycentric)
    ismissing(r.stats) && error("No convergence data")
    fig = Figure( )
    ax = Axis(
        fig[1,1], xlabel="degree", ylabel="max error", yscale=log10
        )
    stats = r.stats
    deg = length.(stats.nodes) .- 1
    scatter!(ax, deg, stats.error,
        color=(stats.nbad.==0), colormap=[:red, :darkblue], colorrange=(0, 1), markersize=10)
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
    t, τ = RFA.refine(p, r.prenodes, max(30, ceil(Int,600/degree(r))))
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

end # module
