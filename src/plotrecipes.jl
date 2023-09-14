# Plot a convergence curve, showing which steps had bad poles
@recipe function foo(s::ConvergenceStats; goodcolor=:darkblue, badcolor=:red)
    seriestype := :path
    linecolor --> :black
    markershape --> :circle
    markerstrokewidth --> 0
    yscale --> :log10
    ylabel --> "max error"
    xlabel --> "degree"
    label --> nothing
    markercolor := [ n > 0 ? badcolor : goodcolor for n in s.nbad ]
    s.error
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


# function convergenceplot(stats::ConvergenceStats)
#     fig = Figure( resolution=(800, 800) )
#     ax = Axis(
#         fig[1,1], 
#         xlabel="degree", ylabel="error", yscale=log10
#         )
#     data = Point.(eachindex(stats.error), stats.error)
#     scatter!(ax, data, color=(stats.nbad.==0), colormap=[:red, :darkblue], colorrange=(0, 1), markersize=10)
#     scatter!(ax, data[stats.bestidx], color=RGBAf(1,1,1,0), strokecolor=:yellow, strokewidth=3, markersize=16 )

#     x = sort(last(stats.nodes))
#     ax = Axis(
#         fig[2,1], 
#         xlabel="x", ylabel="grid spacing", yscale=log10
#         )
#     lines!(ax, x[1:end-1], diff(x))
#     fig
# end

# function convergenceplot(f::Function, r::Barycentric, stats::ConvergenceStats)
#     fig = Figure( resolution=(800, 800) )
#     ax = Axis(
#         fig[1,1], 
#         xlabel="degree", ylabel="error", yscale=log10
#         )
#     data = Point.(eachindex(stats.error), stats.error)
#     scatter!(ax, data, color=(stats.nbad.==0), colormap=[:red, :darkblue], colorrange=(0, 1), markersize=10)
#     scatter!(ax, data[stats.bestidx], color=RGBAf(1,1,1,0), strokecolor=:yellow, strokewidth=3, markersize=16 )

#     x = stats.nodes[stats.bestidx]
#     xx = sort(refine(x, 30))

#     ax = Axis(
#         fig[2,1], 
#         xlabel="x", ylabel="error"
#         )
#     lines!(ax, xx, f.(xx) - r.(xx))
#     fig
# end
