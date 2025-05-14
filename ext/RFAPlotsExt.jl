module RFAPlotsExt

using RationalFunctionApproximation, ComplexRegions, Plots
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

function RFA.errorplot(r::RFA.Approximation; use_abs=false)
    fig = plot(xlabel="boundary parameter", ylabel="error", legend=false)
    # try to get enough points for a smooth result
    N = ceil(Int, 1000 / length(nodes(r)))
    t, τ, err = check(r, refinement=N, quiet=true, prenodes=true)
    if use_abs
        plot!(t, abs.(err))
        ylabel!("| error |")
    elseif isreal(err)
        plot!(t, err)
        ylabel!("error")
    else
        plot!(t, real.(err))
        plot!(t, imag.(err))
        ylabel!("Re, Im error")
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
    fig = plot(xlabel="Re(z)", ylabel="Im(z)", aspect_ratio=1, legend=false)
    z, _ = check(r, quiet=true)
    z = [z; z[1]]
    plot!(real(z), imag(z))
    zp = iszero(idx) ? RFA.poles(r) : RFA.poles(rewind(r, idx))
    color = [r.allowed(z) ? :black : :red for z in zp]
    scatter!(real(zp), imag(zp); color, markersize=4)
    xbox, ybox = axisbox(z)
    xlims!(xbox...)
    ylims!(ybox...)
    return fig
end

end # module
