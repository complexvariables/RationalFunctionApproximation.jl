using RationalFunctionApproximation, ComplexRegions, CairoMakie, LaTeXStrings
using Logging, Printf

const RFA = RationalFunctionApproximation
isdefined(RFA, :QuadraticThiele) || error(
    "These experiments require the local QTCF extension. Run this file with --project=test from this checkout.")

const ROOT = @__DIR__

const OUTPUT = joinpath(ROOT, "qtcf_output")
const FIGURES = joinpath(OUTPUT, "figures")
const NUMBERS = joinpath(OUTPUT, "numbers")
mkpath.((FIGURES, NUMBERS))

const T = Float64
const DEFAULT_TOL = 1e-12
const AAA_COLOR = RGBf(0 / 255, 114 / 255, 178 / 255)
const TCF_COLOR = RGBf(230 / 255, 159 / 255, 0 / 255)
const QTCF_COLOR = RGBf(0 / 255, 158 / 255, 115 / 255)

global_logger(ConsoleLogger(stderr, Logging.Error))
set_eval_method(OneDiv())
set_weight_method(OneDiv())

function save_figure(name, figure)
    save(joinpath(FIGURES, name * ".pdf"), figure)
    save(joinpath(FIGURES, name * ".png"), figure; px_per_unit=2)
end

println("RationalFunctionApproximation ", pkgversion(RFA))
println("Package source: ", pathof(RFA))
println("Running all Section 6 experiments with full paper settings.")
println("Output: ", OUTPUT)

function history_data(result, f, points)
    target = f.(points)
    approximants = getproperty.(result.history, :interpolant)
    denominator_degrees = [degrees(approximant)[2] for approximant in approximants]
    errors = Float64[]
    for approximant in approximants
        pointwise = abs.(target .- approximant.(points))
        push!(errors, all(isfinite, pointwise) ? maximum(pointwise) : Inf)
    end
    candidates = findall(isfinite, errors)
    isempty(candidates) && error("No finite interval approximation was found.")
    best = candidates[argmin(errors[candidates])]
    (; degrees=denominator_degrees, errors, best)
end

fine = T(2) .^ (-100:0.1:-1)
abs_check = sort!(unique(vcat(
    range(-one(T), one(T), 10001),
    fine,
    -fine,
    fine .- one(T),
)))
abs_f = x -> abs(real(x))
abs_domain = Segment{T}(-1, 1)
abs_symmetric = RFA.Symmetric(abs_domain; reflection=x -> -conj(x))

tcf_max_iter = 241
tcf_parameters = (
    tol=100eps(T), allowed=true,
    max_iter=tcf_max_iter,
    refinement=3,
    stagnation=tcf_max_iter,
)
qtcf_parameters = (
    tol=100eps(T), allowed=true,
    max_iter=150,
    refinement=4,
    initial_refinement=28,
    stagnation=170,
)

println("\nInterval case: |x| on [-1,1]")
set_weight_method(OneDiv())
abs_tcf = approximate(abs_f, abs_domain, TCF(); tcf_parameters...)
set_weight_method(OneDiv())
abs_qtcf_symmetric = approximate(
    abs_f, abs_symmetric, QuadraticThiele(); qtcf_parameters...)

imag_f = z -> abs(z)
imag_domain = Segment(-1im, 1im)
imag_symmetric = RFA.Symmetric(imag_domain; reflection=conj)
imag_check = im .* abs_check

println("Interval case: |z| on [-i,i]")
set_weight_method(OneDiv())
imag_tcf = approximate(imag_f, imag_domain, TCF(); tcf_parameters...)
set_weight_method(OneDiv())
imag_qtcf_symmetric = approximate(
    imag_f, imag_symmetric, QuadraticThiele(); qtcf_parameters...)

interval_results = (
    real=(
        title=raw"$|x|$ on $[-1,1]$",
        f=abs_f,
        points=abs_check,
        runs=(
            (label="TCF", result=abs_tcf),
            (label="symmetric QTCF", result=abs_qtcf_symmetric),
        ),
    ),
    imaginary=(
        title=raw"$|z|$ on $[-i,i]$",
        f=imag_f,
        points=imag_check,
        runs=(
            (label="TCF", result=imag_tcf),
            (label="symmetric QTCF", result=imag_qtcf_symmetric),
        ),
    ),
)

function interval_outputs()
    figure = Figure(size=(720, 340), fontsize=11)
    panels = (interval_results.real, interval_results.imaginary)
    colors = (TCF_COLOR, QTCF_COLOR)
    labels = ("TCF", "symmetric QTCF")

    open(joinpath(NUMBERS, "interval_convergence.csv"), "w") do io
        println(io, "domain,method,iteration,denominator_degree,max_error,best")
        for (column, panel) in enumerate(panels)
            axis = Axis(figure[1, column]; yscale=log10,
                title=LaTeXString(panel.title), xlabel="denominator degree",
                ylabel="max error", titlesize=14)
            for (run, color, label) in zip(panel.runs, colors, labels)
                data = history_data(run.result, panel.f, panel.points)
                errors = max.(data.errors, floatmin(T))
                scatterlines!(axis, data.degrees, errors;
                    color, label, linewidth=1.0, markersize=2.8, alpha=0.65)
                scatter!(axis, [data.degrees[data.best]], [errors[data.best]];
                    color=:transparent, strokecolor=:black,
                    strokewidth=1.8, markersize=10)
                domain_name = column == 1 ? "real" : "imaginary"
                for j in eachindex(errors)
                    println(io, join((domain_name, run.label, j,
                        data.degrees[j], data.errors[j], j == data.best), ','))
                end
            end
            text!(axis, 0.96, 0.94; text=column == 1 ? "(a)" : "(b)",
                space=:relative, align=(:right, :top), fontsize=14)
        end
    end

    legend_elements = [
        [LineElement(color=color, linewidth=1.3),
            MarkerElement(color=color, marker=:circle, markersize=4)]
        for color in colors
    ]
    Legend(figure[2, 1:2], legend_elements, collect(labels);
        orientation=:horizontal, framevisible=true, padding=(4, 4, 2, 2),
        labelsize=10)
    rowgap!(figure.layout, 2)
    colgap!(figure.layout, 12)
    save_figure("Figure_6_1_interval_convergence", figure)
end

interval_outputs()

const SCHWARZ_AAA_PARAMETERS = (
    tol=DEFAULT_TOL,
    allowed=true,
    max_iter=150,
    refinement=4,
    stagnation=10,
)
const SCHWARZ_TCF_PARAMETERS = (
    tol=DEFAULT_TOL,
    allowed=true,
    max_iter=500,
    refinement=3,
    stagnation=10,
)
const SCHWARZ_QTCF_PARAMETERS = (
    tol=DEFAULT_TOL,
    allowed=true,
    max_iter=100,
    refinement=3,
    initial_refinement=20,
    stagnation=10,
)

raw(result) = get_function(result)
curve_points(domain, n=800) = [point(domain, t) for t in range(0, length(domain); length=n)]
grid_points(domain, n=1200) = RFA.isclosed(domain) ?
    [point(domain, t) for t in range(0, length(domain); length=n + 1)[1:end-1]] :
    [point(domain, t) for t in range(0, length(domain); length=n)]

function finite_max_error(result, f, points)
    result === nothing && return NaN
    values = abs.(f.(points) .- result.(points))
    any(!isfinite, values) && return Inf
    maximum(values; init=0.0)
end

half_turn(c) = z -> 2c - z
line_reflection(c, factor) = z -> begin
    w = z - c
    c + factor * complex(real(w), -imag(w))
end

upper_ellipse_point(t) = begin
    theta = pi * (1 - t)
    cos(theta) + 0.5im * sin(theta)
end
upper_ellipse_deriv(t) = begin
    theta = pi * (1 - t)
    pi * sin(theta) - 0.5im * pi * cos(theta)
end
paper_upper_ellipse = Curve(upper_ellipse_point, upper_ellipse_deriv, 0, 1)

paper_lobed_full = ClosedCurve(t -> (1 + 0.2sin(5t)) * cis(t), 0, 2pi)

paper_smooth_s_point(t) = -0.75sin(2pi * t) + im * (2 - 4t)
paper_smooth_s_deriv(t) = -1.5pi * cos(2pi * t) - 4im
paper_smooth_s_curve = Curve(paper_smooth_s_point, paper_smooth_s_deriv, 0, 1)

function paper_semicircle_s_point(t)
    if t <= 0.5
        s = 2t
        return im + cis(pi / 2 + pi * s)
    end
    s = 2t - 1
    -im + cis(pi / 2 - pi * s)
end

paper_semicircle_s_deriv(t) =
    (paper_semicircle_s_point(t + 1e-6) - paper_semicircle_s_point(t - 1e-6)) / 2e-6
paper_semicircle_s_curve =
    Curve(paper_semicircle_s_point, paper_semicircle_s_deriv, 0, 1)

signedpow(u, p) = abs(u) < 1e-14 ? 0.0 : sign(u) * abs(u)^p
function super6_point(t)
    signedpow(cos(t), 1 / 3) + im * signedpow(sin(t), 1 / 3)
end
super6_deriv(t) = (super6_point(t + 1e-6) - super6_point(t - 1e-6)) / 2e-6
paper_superellipse6 = ClosedCurve(super6_point, super6_deriv, 0, 2pi)

const inlet_a0 = 0.1
const inlet_A = [-0.1, -0.7, -0.1]
const inlet_B = [-0.8, -0.4, -0.2]

paper_inlet_point(t) = inlet_a0 + sum(
    inlet_A[k] * cos(k * t) + im * inlet_B[k] * sin(k * t)
    for k in eachindex(inlet_A))
paper_inlet_deriv(t) = sum(
    -k * inlet_A[k] * sin(k * t) + im * k * inlet_B[k] * cos(k * t)
    for k in eachindex(inlet_A))
paper_inlet = ClosedCurve(paper_inlet_point, paper_inlet_deriv, 0, 2pi)

origin = 0.0 + 0.0im

function qtcf_case(slug, name, domain, reflection, formula;
    audit_n=1800, curve_n=700)
    (; slug, name, domain, reflection, formula, audit_n, curve_n)
end

qtcf_cases = [
    qtcf_case("ellipse_1_025", "ellipse (1, 1/4) Schwarz", Shapes.ellipse(1, 0.25),
        half_turn(origin), raw"$z(t)=\cos t+\frac{i}{4}\sin t$"),
    qtcf_case("ellipse_1_05", "ellipse (1, 1/2) Schwarz", Shapes.ellipse(1, 0.5),
        half_turn(origin), raw"$z(t)=\cos t+\frac{i}{2}\sin t$"),
    qtcf_case("upper_ellipse_1_05", "upper ellipse (1, 1/2) Schwarz",
        paper_upper_ellipse, line_reflection(origin, -1.0 + 0.0im),
        raw"$z(t)=\cos(\pi(1-t))+\frac{i}{2}\sin(\pi(1-t))$"),
    qtcf_case("five_lobed", "five-lobed Schwarz", paper_lobed_full,
        line_reflection(origin, -1.0 + 0.0im),
        raw"$z(t)=(1+0.2\sin 5t)e^{it}$"),
    qtcf_case("smooth_s", "smooth letter S Schwarz", paper_smooth_s_curve,
        half_turn(origin), raw"$z(t)=-0.75\sin(2\pi t)+i(2-4t)$"),
    qtcf_case("two_semicircle_s", "two-semicircle S Schwarz", paper_semicircle_s_curve,
        half_turn(origin),
        raw"$z(t)=i+e^{i(\pi/2+2\pi t)}\ (0\leq t\leq1/2),\quad z(t)=-i+e^{i(3\pi/2-2\pi t)}\ (1/2<t\leq1)$"),
    qtcf_case("superellipse_6", "superellipse x^6+y^6=1 Schwarz", paper_superellipse6,
        half_turn(origin),
        raw"$z(t)=\mathrm{sgn}(\cos t)|\cos t|^{1/3}+i\mathrm{sgn}(\sin t)|\sin t|^{1/3}$"),
    qtcf_case("analytic_inlet", "analytic right inlet Schwarz", paper_inlet,
        line_reflection(origin, 1.0 + 0.0im),
        raw"$z(t)=0.1-0.1\cos t-0.7\cos2t-0.1\cos3t-i(0.8\sin t+0.4\sin2t+0.2\sin3t)$"),
]

function run_method(case, method)
    f = conj
    selector = method === :AAA ? AAA() : method === :TCF ? TCF() : QuadraticThiele()
    start = time_ns()
    result = if method === :QTCF
        symmetric = RFA.Symmetric(case.domain; reflection=case.reflection)
        approximate(f, symmetric, selector; SCHWARZ_QTCF_PARAMETERS...)
    elseif method === :AAA
        approximate(f, case.domain, selector; SCHWARZ_AAA_PARAMETERS...)
    else
        approximate(f, case.domain, selector; SCHWARZ_TCF_PARAMETERS...)
    end
    elapsed = (time_ns() - start) / 1e9
    rational = raw(result)
    pole, residue = residues(rational)
    audit = grid_points(case.domain, case.audit_n)
    error = finite_max_error(result, f, audit)
    rational_degrees = degrees(rational)
    (; result, rational, pole, residue, error, elapsed,
        rational_degree=rational_degrees[2], rational_degrees)
end

function run_case(case)
    println("\nSchwarz case: ", case.name)
    rows = Dict{Symbol,Any}()
    for method in (:AAA, :TCF, :QTCF)
        print("  continuum ", method, " ... ")
        try
            row = run_method(case, method)
            rows[method] = merge(row, (failed=false, exception=nothing))
            @printf("error %.3e, degree %s, %.3f s\n",
                row.error, string(row.rational_degrees), row.elapsed)
        catch exception
            rows[method] = (
                result=nothing, rational=nothing,
                pole=ComplexF64[], residue=ComplexF64[],
                error=NaN, elapsed=NaN, rational_degree=-1,
                rational_degrees=(-1, -1), failed=true, exception=exception)
            println("failed: ", exception)
        end
    end
    (; case, curve=curve_points(case.domain, case.curve_n), rows)
end

schwarz_results = run_case.(qtcf_cases)

const CASE_LABELS = Dict(
    "ellipse_1_05" => "Ellipse (1, 1/2)",
    "ellipse_1_025" => "Ellipse (1, 1/4)",
    "upper_ellipse_1_05" => "Upper ellipse (1, 1/2)",
    "five_lobed" => "Five-lobed curve",
    "smooth_s" => "Smooth S",
    "two_semicircle_s" => "Two-semicircle S",
    "superellipse_6" => "Superellipse",
    "analytic_inlet" => "Analytic inlet",
)

const CONTINUUM_LIMITS = (-2.25, 2.25, -2.25, 2.25)
const CONTINUUM_TICKS = -2:1:2

function write_continuum_summary()
    open(joinpath(NUMBERS, "schwarz_summary.csv"), "w") do io
        println(io, "case,method,max_error,degree,numerator_degree,denominator_degree,time_seconds,failed")
        for result in schwarz_results, method in (:AAA, :TCF, :QTCF)
            row = result.rows[method]
            pdeg, qdeg = row.rational_degrees
            println(io, join((result.case.slug, method, row.error,
                row.rational_degree, pdeg, qdeg, row.elapsed, row.failed), ','))
        end
    end
end

function in_window(poles, limits)
    xmin, xmax, ymin, ymax = limits
    isfinite.(poles) .&
        (xmin .<= real.(poles)) .& (real.(poles) .<= xmax) .&
        (ymin .<= imag.(poles)) .& (imag.(poles) .<= ymax)
end

const CONTINUUM_COLORRANGE = let
    logs = Float64[]
    for result in schwarz_results, method in (:AAA, :TCF, :QTCF)
        row = result.rows[method]
        keep = in_window(row.pole, CONTINUUM_LIMITS) .& isfinite.(row.residue)
        append!(logs, log10.(abs.(row.residue[keep]) .+ eps(T)))
    end
    sort!(logs)
    if isempty(logs)
        lo, hi = -4, 0
    else
        n = length(logs)
        lo = floor(Int, logs[clamp(ceil(Int, 0.10n), 1, n)])
        hi = ceil(Int, logs[clamp(ceil(Int, 0.99n), 1, n)])
        lo == hi && (lo -= 1)
    end
    (Float64(lo), Float64(hi))
end

function residue_power_ticks(lo, hi)
    lower = floor(Int, lo)
    upper = ceil(Int, hi)
    step = max(1, cld(abs(lower), 4))
    exponents = collect(lower:step:-1)
    lo <= 0 <= hi && push!(exponents, 0)
    upper > 0 && push!(exponents, upper)
    unique!(sort!(exponents))
    Float64.(exponents), [LaTeXString("10^{$exponent}") for exponent in exponents]
end

function write_continuum_pole_data(result)
    case = result.case
    open(joinpath(NUMBERS, case.slug * "_poles_residues.csv"), "w") do io
        println(io, "method,index,pole_real,pole_imag,residue_real,residue_imag,log10_abs_residue,in_plot_window")
        for method in (:AAA, :TCF, :QTCF)
            row = result.rows[method]
            count = min(length(row.pole), length(row.residue))
            inside = in_window(row.pole, CONTINUUM_LIMITS)
            for k in 1:count
                residue_log = log10(abs(row.residue[k]) + eps(T))
                println(io, join((method, k,
                    real(row.pole[k]), imag(row.pole[k]),
                    real(row.residue[k]), imag(row.residue[k]),
                    residue_log, inside[k]), ','))
            end
        end
    end
end

function continuum_pole_figure(results, number)
    methods = (:AAA, :TCF, :QTCF)
    lo, hi = CONTINUUM_COLORRANGE

    nrows = length(results)
    figure = Figure(size=(720, 760), fontsize=8.5, figure_padding=2)
    for (j, method) in enumerate(methods)
        Label(figure[0, j], String(method); fontsize=13, font=:bold)
    end
    last_scatter = nothing
    for (i, result) in enumerate(results)
        case = result.case
        plot_row = 2i - 1
        for (j, method) in enumerate(methods)
            row = result.rows[method]
            axis = Axis(figure[plot_row, j]; aspect=DataAspect(), limits=CONTINUUM_LIMITS,
                title=CASE_LABELS[case.slug], titlesize=9,
                xlabel="Re z", ylabel="Im z", xlabelsize=8, ylabelsize=8,
                xticks=CONTINUUM_TICKS, yticks=CONTINUUM_TICKS,
                xticklabelsize=7, yticklabelsize=7)
            lines!(axis, real.(result.curve), imag.(result.curve);
                color=AAA_COLOR, linewidth=1.15)
            keep = in_window(row.pole, CONTINUUM_LIMITS) .& isfinite.(row.residue)
            if any(keep)
                last_scatter = scatter!(axis,
                    real.(row.pole[keep]), imag.(row.pole[keep]);
                    color=log10.(abs.(row.residue[keep]) .+ eps(T)),
                    colormap=Reverse(:thermal), colorrange=(lo, hi), markersize=5.0)
            end
        end
        Label(figure[2i, 1:3], LaTeXString(case.formula); fontsize=9)
    end
    tick_positions, tick_labels = residue_power_ticks(lo, hi)
    last_scatter !== nothing && Colorbar(figure[1:2nrows, 5], last_scatter;
        label="residue", ticks=(tick_positions, tick_labels),
        labelsize=9, ticklabelsize=7)
    for column in 1:3
        colsize!(figure.layout, column, Fixed(175))
    end
    colsize!(figure.layout, 4, Fixed(24))
    for row in 1:nrows
        rowsize!(figure.layout, 2row - 1, Fixed(175))
        rowsize!(figure.layout, 2row, Fixed(17))
    end
    rowgap!(figure.layout, 2)
    for row in 1:nrows-1
        rowgap!(figure.layout, 2row, 10)
    end
    colgap!(figure.layout, 16)
    colgap!(figure.layout, 3, 0)
    colgap!(figure.layout, 4, 0)
    resize_to_layout!(figure)

    name = "Figure_6_$(number)_continuum_poles_$(number - 1)"
    save_figure(name, figure)
    name
end

write_continuum_summary()
write_continuum_pole_data.(schwarz_results)
groups = [schwarz_results[i:min(i + 3, end)] for i in 1:4:length(schwarz_results)]
for (i, group) in enumerate(groups)
    continuum_pole_figure(group, i + 1)
end

failed_runs = [(result.case.slug, method) for result in schwarz_results
    for method in (:AAA, :TCF, :QTCF) if result.rows[method].failed]
isempty(failed_runs) || error("Some Schwarz experiments failed; inspect numbers/schwarz_summary.csv: $failed_runs")

println("Completed all 28 fits for Figures 6.1–6.3.")
println("Saved 3 figures as PDF and PNG, plus 10 CSV data files, to ", OUTPUT)
