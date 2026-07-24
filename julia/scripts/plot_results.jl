#!/usr/bin/env julia

using CSV
using JSON3
using Plots

function options(args)
    parsed = Dict{String,String}()
    index = 1
    while index <= length(args)
        startswith(args[index], "--") || error("unexpected argument: $(args[index])")
        if args[index] == "--double-dimer"
            parsed["double-dimer"] = "true"
            index += 1
            continue
        end
        index < length(args) || error("missing value for $(args[index])")
        parsed[args[index][3:end]] = args[index + 1]
        index += 2
    end
    return parsed
end

function model_label(distribution, params)
    parsed = JSON3.read(String(params))
    isempty(parsed) && return String(distribution)
    values = join(("$key=$(round(Float64(value); digits=3))" for (key, value) in pairs(parsed)), ", ")
    return "$(replace(String(distribution), "_" => " ")) ($values)"
end

function safe_name(distribution, params)
    label = lowercase(model_label(distribution, params))
    return strip(replace(label, r"[^a-z0-9]+" => "_"), '_')
end

function scaling_panel(summary_rows, model_rows, title; double_dimer=false)
    ordered = sort(summary_rows; by=row -> row.L)
    xs = Float64.(getproperty.(ordered, :log_L))
    ys = Float64.(getproperty.(ordered, :annealed_variance))
    errors = 1.96 .* Float64.(getproperty.(ordered, :annealed_variance_se))
    panel = scatter(xs, ys; yerror=errors, markercolor=:black, markerstrokewidth=0,
                    markersize=3.5, label="simulation", title, xlabel="log L",
                    ylabel=double_dimer ? "Var(W₁ - W₂)" : "Var(W_L)",
                    legend=:topleft, legendfontsize=7,
                    titlefontsize=9, guidefontsize=8, tickfontsize=7)
    grid = range(minimum(xs), maximum(xs); length=200)
    for row in model_rows
        if row.model == "a_plus_b_log_L"
            plot!(panel, grid, row.intercept .+ row.coefficient .* grid;
                  color=:steelblue, linewidth=2, label="a + b log L")
        elseif row.model == "a_plus_b_log_L_squared"
            plot!(panel, grid, row.intercept .+ row.coefficient .* grid .^ 2;
                  color=:darkorange, linewidth=2, linestyle=:dash,
                  label="a + b (log L)^2")
        end
    end
    return panel
end

function main(args=ARGS)
    parsed = options(args)
    analysis_dir = get(parsed, "analysis-dir", "analysis_strict_annealed")
    output_dir = get(parsed, "output-dir", "../reports")
    double_dimer = haskey(parsed, "double-dimer")
    figures_dir = joinpath(output_dir, "figures")
    details_dir = joinpath(figures_dir, "by_distribution")
    mkpath(details_dir)

    summary = NamedTuple.(CSV.File(joinpath(analysis_dir, "summary.csv")))
    fits = NamedTuple.(CSV.File(joinpath(analysis_dir, "loglog_fits.csv")))
    models = NamedTuple.(CSV.File(joinpath(analysis_dir, "scaling_model_comparison.csv")))
    keys = sort!(unique((String(row.distribution), String(row.distribution_params))
                        for row in summary))

    panels = Any[]
    for key in keys
        summary_rows = [row for row in summary
                        if String(row.distribution) == key[1] && String(row.distribution_params) == key[2]]
        model_rows = [row for row in models
                      if String(row.distribution) == key[1] &&
                         String(row.distribution_params) == key[2] &&
                         row.variance_kind == "annealed"]
        title = model_label(key...)
        panel = scaling_panel(summary_rows, model_rows, title; double_dimer)
        push!(panels, panel)
        detail = scaling_panel(summary_rows, model_rows, title; double_dimer)
        savefig(detail, joinpath(details_dir, safe_name(key...) * ".png"))
        savefig(detail, joinpath(details_dir, safe_name(key...) * ".pdf"))
    end
    overview_title = double_dimer ? "Double-dimer winding-difference variance" :
                                   "Strict-annealed winding variance"
    overview = plot(panels...; layout=(4, 4), size=(1600, 1400),
                    plot_title=overview_title)
    savefig(overview, joinpath(figures_dir, "annealed_scaling_all_distributions.png"))
    savefig(overview, joinpath(figures_dir, "annealed_scaling_all_distributions.pdf"))

    ordered_fits = sort(fits; by=row -> row.p)
    labels = model_label.(getproperty.(ordered_fits, :distribution),
                          getproperty.(ordered_fits, :distribution_params))
    ps = Float64.(getproperty.(ordered_fits, :p))
    lows = Float64.(getproperty.(ordered_fits, :p_bootstrap_ci_low))
    highs = Float64.(getproperty.(ordered_fits, :p_bootstrap_ci_high))
    forest = scatter(ps, eachindex(ps); xerror=(ps .- lows, highs .- ps),
                     yticks=(eachindex(ps), labels), markercolor=:black,
                     markerstrokewidth=0, label="95% bootstrap CI", xlabel="p",
                     title=double_dimer ? "Var(W₁ - W₂) = C (log L)^p" :
                                          "Var(W_L) = C (log L)^p",
                     size=(1000, 800),
                     left_margin=8Plots.mm, legend=:bottomright)
    vline!(forest, [1.0]; color=:seagreen, linewidth=2, label="p = 1")
    vline!(forest, [2.0]; color=:firebrick, linewidth=2, linestyle=:dash, label="p = 2")
    savefig(forest, joinpath(figures_dir, "scaling_exponent_forest.png"))
    savefig(forest, joinpath(figures_dir, "scaling_exponent_forest.pdf"))

    bic_labels = String[]
    deltas = Float64[]
    for key in keys
        rows = [row for row in models if String(row.distribution) == key[1] &&
                String(row.distribution_params) == key[2] && row.variance_kind == "annealed"]
        log_row = only(row for row in rows if row.model == "a_plus_b_log_L")
        squared_row = only(row for row in rows if row.model == "a_plus_b_log_L_squared")
        push!(bic_labels, model_label(key...))
        push!(deltas, squared_row.bic - log_row.bic)
    end
    ordering = sortperm(deltas)
    bic = bar(deltas[ordering], orientation=:h,
              yticks=(eachindex(ordering), bic_labels[ordering]),
              color=[value >= 0 ? :seagreen : :firebrick for value in deltas[ordering]],
              label="", xlabel="BIC(log^2) - BIC(log)",
              title="Positive values favor ordinary logarithmic growth",
              size=(1000, 800), left_margin=8Plots.mm)
    vline!(bic, [0.0]; color=:black, linewidth=1, label="")
    savefig(bic, joinpath(figures_dir, "bic_model_comparison.png"))
    savefig(bic, joinpath(figures_dir, "bic_model_comparison.pdf"))

    for name in ("summary.csv", "loglog_fits.csv", "scaling_model_comparison.csv",
                 "pointwise_ratios.csv", "local_effective_exponents.csv")
        cp(joinpath(analysis_dir, name), joinpath(output_dir, name); force=true)
    end
    pair_path = joinpath(analysis_dir, "double_dimer_pairs.csv")
    isfile(pair_path) && cp(pair_path, joinpath(output_dir, "double_dimer_pairs.csv"); force=true)
    println("Wrote Julia reports to $output_dir")
end

main()
