module RFAPlotsExt

using RationalFunctionApproximation, Plots
RFA = RationalFunctionApproximation

# Plot a convergence curve, showing which steps had bad poles
function RFA.convergenceplot(r::RFA.Approximation)
    deg, err, _, allowed, best = get_history(r)
    fig = plot(xlabel="degree", ylabel="max relative error", yscale=:log10, legend=false)
    seriescolor = [ all(b) ? :darkblue : :red for b in allowed]
    scatter!([deg[best]], [err[best]],
        markercolor=RGBA(1,1,1,0), markerstrokecolor=:gold, msw=5, markersize=7 )
    scatter!(deg, err; seriescolor, msw=0, markersize=5)
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
