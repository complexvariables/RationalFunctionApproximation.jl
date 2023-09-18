module RFAPlotsExt

using RationalFunctionApproximation, Plots
RFA = RationalFunctionApproximation

# Plot a convergence curve, showing which steps had bad poles
RFA.convergenceplot(r::RFA.Approximation) = RFA.convergenceplot(r.fun)
function RFA.convergenceplot(r::Barycentric)
    ismissing(r.stats) && error("No convergence data")
    fig = plot(xlabel="degree", ylabel="max error", yscale=:log10, legend=false)
    stats = r.stats
    deg = length.(stats.nodes) .- 1
    seriescolor = [ n > 0 ? :red : :darkblue for n in stats.nbad]
    scatter!([deg[stats.bestidx]], [stats.error[stats.bestidx]],
    markercolor=RGBA(1,1,1,0), markerstrokecolor=:gold, msw=5, markersize=7 )
    scatter!(deg, stats.error; seriescolor, msw=0, markersize=5)
    return fig
end

# Plot the domain and the poles of the approximation
@userplot PolePlot
@recipe function foo(PP::PolePlot)
    @assert length(PP.args) == 1 && (typeof(PP.args[1]) <: Approximation)
        "Must be given an Approximation object."

    r = PP.args[1]

    @series begin
        label --> nothing
        aspect_ratio --> 1
        r.domain
    end

    @series begin
        zp = poles(r.fun)
        seriestype := :scatter
        color --> :red
        markerstrokewidth --> 0
        label --> nothing
        aspect_ratio --> 1
        real(zp), imag(zp)
    end
end


end # module
