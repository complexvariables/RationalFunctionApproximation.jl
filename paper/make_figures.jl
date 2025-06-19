using RationalFunctionApproximation
f = z -> log(1 + 1im + 5im*z)
r = approximate(f, unit_interval)

using CairoMakie, LaTeXStrings
poleplot(r)
save("paper/poleplot1.pdf", current_figure())

##
f = z -> sqrt(z + 1e-6im)
x = range(-1, 1, 1001)
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel=L"x", ylabel=L"\mathrm{Re}(f(x) - r(x))")
r = approximate(f, x)
lines!(ax1, -0.001..0.001, x->real(f(x) - r(x)))

r = approximate(f, unit_interval)
ax2 = Axis(fig[2, 1], xlabel=L"x", ylabel=L"\mathrm{Re}(f(x) - r(x))")
lines!(ax2, -0.001..0.001, x->real(f(x) - r(x)))
save("paper/continuum.pdf", current_figure())

##
@elapsed r = approximate(abs, unit_interval)
# 0.229 sec
deg, err, _, allowed, best = get_history(r);
fig = Figure( )
color = [all(b) ? :darkblue : :red for b in allowed]
ax = Axis(fig[1,1], xlabel="degree", ylabel="relative max error", yscale=log10)
scatter!(ax, deg, err; color, colorrange=(0, 1), markersize=8)
scatter!(ax, deg[best], err[best],
    color=RGBAf(1,1,1,0), strokecolor=:gold, strokewidth=3, markersize=13 )

##
using DoubleFloats, ComplexRegions
uiquad = Segment{Double64}(-1, 1)
@elapsed rq = approximate(abs, uiquad; allowed=true, max_iter=200)
# 30.46 sec

##
deg = [degree(h.interpolant) for h in rq.history]
err = [h.error for h in rq.history]
ax = Axis(fig[2,1], xlabel="degree", xscale=sqrt, ylabel="relative max error", yscale=log10)
scatter!(ax, deg[2:end], err[2:end]; color=:black, markersize=8)
save("paper/absconverge.pdf", current_figure())
fig

##
fig = Figure()
x = sort(filter(>(0), nodes(rq)))
ax = Axis(fig[1,1], xlabel="node index", xscale=sqrt, ylabel="distance to singularity", yscale=log10)
scatter!(ax, x; markersize=8)
save("paper/absnodes.pdf", current_figure())

##
#@elapsed rt = approximate(abs, uiquad; allowed=true, method=Thiele, max_iter=400)
# 1.69 seconds
deg = [degree(h.interpolant) for h in rt.history]
err = [h.error for h in rt.history]
fig = Figure()
ax = Axis(fig[1,1], xlabel="degree", xscale=sqrt, ylabel="relative max error", yscale=log10)
scatter!(ax, deg[3:end], err[3:end]; color=:black, markersize=8)
save("paper/absthiele.pdf", current_figure())


##
f = z -> abs(z - 0.5 + 0.05im)
r = approximate(f, unit_interval, max_iter=20)

fig = Figure()
ax1 = Axis(fig[1,1], xlabel=L"x", ylabel=L"f(x) - r(x)")
lines!(ax1, check(r)...)

r_inf = minimax(r, 20);
ax2 = Axis(fig[2,1], xlabel=L"x", ylabel=L"f(x) - r(x)")
lines!(ax2, check(r_inf)...)
save("paper/minimax.pdf", fig)

##
using SpecialFunctions: zeta
f = z -> 1/ zeta(11z)
r = approximate(f, unit_circle)
poleplot(r)
save("paper/zetapoles.pdf", current_figure())

##
f = z -> coth(1/z^3)
r = approximate(f, exterior(Shapes.squircle))
poleplot(r)
save("paper/squircle.pdf", current_figure())
