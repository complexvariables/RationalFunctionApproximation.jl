function foo(curve::RFA.ComplexCurveOrPath, s::AbstractVector; refinement=0, maxpoints=length(s))
    F = real_type(curve)
    param = Matrix{F}(undef, maxpoints, refinement+1)
    n = length(s)
    for i in 1:n-1
        for j in 0:refinement
            param[i, j+1] = s[i]
        end
    end
    param[n, 1] = s[end]
    points = Matrix{complex(F)}(undef, maxpoints, refinement+1)
    for i in 1:n-1, j in 1:refinement+1
        points[i, j] = curve(param[i, j])
    end
    # points[1:n-1, :] .= point(curve, view(param,1:n-1, :))
    points[n, 1] = curve(s[end])
end


function bar(c)
    k = 1
    p = []
    while k < maximum(c.next)
        push!(p,c.points[k,1])
        k = c.next[k];
    end
    return p
end
