c = ClosedCurve(t -> cispi(2t) + 0.2cospi(6t) - 0.1sinpi(4t))
@test approximate(tan, c) isa RFA.Approximation
@test approximate(tanh, interior(c)) isa RFA.Approximation
@test approximate(exp, exterior(c)) isa RFA.Approximation
