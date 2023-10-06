using RationalFunctionApproximation
f = x -> 1/10^2 / sin(x^2 + 2x/10 + 2/10^2)
r = approximate(f, unit_interval)

##
using CairoMakie
update_theme!(figure_resolution = (500, 300),)

##
lines(-1..1, x -> f(x) - r(x),
    axis = (xlabel="x", ylabel="f(x) - r(x)"),
    figure = (resolution = (600, 400),))
scatter!(nodes(r), 0*nodes(r), color=:black)
save("error_sin.pdf", current_figure())
current_figure()

##
domaincolor(r, [-1.5, 1.5, -0.2, 2.8], abs=true,
figure = (resolution = (600, 400),))
lines!([-1,1],[0,0],color=:white,linewidth=4)
current_axis().aspect = DataAspect()
current_axis().xlabel = "Re(z)"
current_axis().ylabel = "Im(z)"
save("complex.pdf", current_figure())
current_figure()

##
rm = minimax(r)
lines(-1..1, x -> f(x) - rm(x),
    axis = (xlabel="x", ylabel="f(x) - r(x)"),
    figure = (resolution = (600, 400),))

save("error_minimax.pdf", current_figure())

##
g = x -> tanh(5000 * (x - 0.2))
x = -1:0.001:1
r = aaa(x, g.(x))
xx = -1:.0001:1
println("error = $(maximum(abs, @. g(xx) - r(xx) ))")

##
r = approximate(g, unit_interval)
println("error = $(maximum(abs, @. g(xx) - r(xx) ))")

##
r = approximate(z -> tanh(1/z^4), exterior(unit_circle))
domaincolor(r, [-1.5, 1.5, -1.5, 1.5], abs=true,
figure = (resolution = (600, 400),))

save("tanh.pdf", current_figure())
