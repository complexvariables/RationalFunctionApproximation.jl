n = 17;
eta = exp(-1/sqrt(n))
xpos = eta.^(n-1:-1:0)
x = [-reverse(xpos); 0; xpos]
t, h = approximate(abs.(x), x, Barycentric; max_iter=2n+1)
findmax(abs, t.(x) - abs.(x))
