@usingany CairoMakie, Colors
ff(x) = sqrt(x^2-1)
# r = approximate(z->log(z+0.05+0.05im), unit_interval)
r = approximate(ff, Segment{Double64}(101//100+1/20im,3//2+1/20im); method=Thiele, allowed=true, tol=1e-30)
jgreen = RGB(0.22, 0.596, 0.149)
jred = RGB(0.796,0.235,0.2)
jpurple = RGB(0.584, 0.345, 0.698)
jblue = RGB(0.255, 0.412, 0.882);

##
z = last(poles(r), 200)
scatter(reim(z)..., axis=(autolimitaspect=1,))

##
fig = Figure(size=(600,600),backgroundcolor=:transparent)
# ax = Axis(fig[1,1], autolimitaspect=1)
ax = Axis(fig[1,1],backgroundcolor=:transparent)
hidespines!(ax)
hidedecorations!(ax)
zz = first(sort(z, by=real, rev=true), 9)
scatter!(ax, reim(zz)...,
    color=repeat(reverse([jgreen,jred,jpurple]),outer=3),
    markersize=[round(Int,12*(1.44)^j) for j in 0:8] )
limits!(ax, 0.912,1.001,-0.01,0.06)
fig
