#!/usr/bin/env julia

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CSV
using JSON3
using Plots

function options(args)
    parsed = Dict{String,String}()
    index = 1
    while index <= length(args)
        startswith(args[index], "--") || error("unexpected argument: $(args[index])")
        if args[index] in ("--paired", "--double-dimer")
            parsed["paired"] = "true"
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
    name = titlecase(replace(String(distribution), "_" => " "))
    isempty(parsed) && return name
    parameter_labels = Dict("alpha" => "α", "sigma" => "σ", "shape" => "k")
    values = join(("$(get(parameter_labels, String(key), String(key)))=" *
                   "$(round(Float64(value); digits=3))"
                   for (key, value) in pairs(parsed)), ", ")
    return "$name ($values)"
end

function safe_name(distribution, params)
    parsed = JSON3.read(String(params))
    pieces = String[String(distribution)]
    append!(pieces, ("$(key)_$(round(Float64(value); digits=3))"
                     for (key, value) in pairs(parsed)))
    label = lowercase(join(pieces, "_"))
    return strip(replace(label, r"[^a-z0-9]+" => "_"), '_')
end

function scaling_panel(summary_rows, model_rows, title; paired=false)
    ordered = sort(summary_rows; by=row -> row.L)
    xs = Float64.(getproperty.(ordered, :log_L))
    ys = Float64.(getproperty.(ordered, :annealed_variance))
    errors = 1.96 .* Float64.(getproperty.(ordered, :annealed_variance_se))
    panel = scatter(xs, ys; yerror=errors, markercolor=:black, markerstrokewidth=0,
                    markersize=3.5, label="simulation", title, xlabel="log L",
                    ylabel=paired ? "Var(ΔWₗ)" : "Var(Wₗ)",
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
                  label="a + b (log L)²")
        end
    end
    return panel
end

function main(args=ARGS)
    parsed = options(args)
    analysis_dir = get(parsed, "analysis-dir", "analysis_strict_annealed")
    output_dir = get(parsed, "output-dir", "../reports")
    paired = haskey(parsed, "paired")
    figures_dir = joinpath(output_dir, "figures")
    details_dir = joinpath(figures_dir, "by_distribution")
    mkpath(details_dir)

    summary = NamedTuple.(CSV.File(joinpath(analysis_dir, "summary.csv")))
    fits = NamedTuple.(CSV.File(joinpath(analysis_dir, "loglog_fits.csv")))
    models = NamedTuple.(CSV.File(joinpath(analysis_dir, "scaling_model_comparison.csv")))
    keys = sort!(unique((String(row.distribution), String(row.distribution_params))
                        for row in summary))
    temporal = !isempty(summary) && hasproperty(first(summary), :environment_model) &&
               all(row -> row.environment_model == "temporal_iid", summary)

    panels = Any[]
    for key in keys
        summary_rows = [row for row in summary
                        if String(row.distribution) == key[1] && String(row.distribution_params) == key[2]]
        model_rows = [row for row in models
                      if String(row.distribution) == key[1] &&
                         String(row.distribution_params) == key[2] &&
                         row.variance_kind == "annealed"]
        title = model_label(key...)
        panel = scaling_panel(summary_rows, model_rows, title; paired)
        push!(panels, panel)
        detail = scaling_panel(summary_rows, model_rows, title; paired)
        savefig(detail, joinpath(details_dir, safe_name(key...) * ".png"))
        savefig(detail, joinpath(details_dir, safe_name(key...) * ".pdf"))
    end
    overview_title = paired ? "Paired LERW winding-difference variance" :
                     temporal ? "Temporal-i.i.d. annealed LERW winding variance" :
                                "Single-walk LERW winding variance"
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
                     title=paired ? "Paired LERWs: Var(ΔWₗ) = C(log L)ᵖ" :
                           temporal ? "Temporal-i.i.d.: Var(Wₗ) = C(log L)ᵖ" :
                                      "Var(Wₗ) = C(log L)ᵖ",
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
              label="", xlabel="ΔBIC = BIC(log²) − BIC(log)",
              title="Additive scaling model comparison",
              size=(1000, 800), left_margin=8Plots.mm,
              right_margin=5Plots.mm, bottom_margin=5Plots.mm)
    vline!(bic, [0.0]; color=:black, linewidth=1, label="")
    savefig(bic, joinpath(figures_dir, "bic_model_comparison.png"))
    savefig(bic, joinpath(figures_dir, "bic_model_comparison.pdf"))

    if temporal
        frequency_panels = Any[]
        direction_columns = (
            ("N", :north_step_frequency), ("E", :east_step_frequency),
            ("S", :south_step_frequency), ("W", :west_step_frequency),
        )
        for key in keys
            rows = sort([row for row in summary
                         if String(row.distribution) == key[1] &&
                            String(row.distribution_params) == key[2]]; by=row -> row.L)
            panel = plot(; xscale=:log2, xlabel="L", ylabel="raw-step frequency",
                         title=model_label(key...), legend=:best, ylim=(0.245, 0.255))
            for (label, column) in direction_columns
                plot!(panel, getproperty.(rows, :L), Float64.(getproperty.(rows, column));
                      marker=:circle, markersize=3, label)
            end
            hline!(panel, [0.25]; color=:black, linestyle=:dash, label="1/4")
            push!(frequency_panels, panel)
        end
        frequencies = plot(frequency_panels...; layout=(2, 2), size=(1200, 900),
                           plot_title="Temporal-i.i.d. raw direction frequencies")
        savefig(frequencies, joinpath(figures_dir, "temporal_direction_frequencies.png"))
        savefig(frequencies, joinpath(figures_dir, "temporal_direction_frequencies.pdf"))

        comparison_path = joinpath(analysis_dir, "temporal_baseline_comparisons.csv")
        comparisons = NamedTuple.(CSV.File(comparison_path))
        comparison_keys = sort!(unique((String(row.distribution),
                                        String(row.distribution_params))
                                       for row in comparisons))
        ratio_plot = plot(; xscale=:log2, xlabel="L", ylabel="variance ratio",
                          title="Temporal distributions / temporal baseline")
        winding_plot = plot(; xscale=:log2, xlabel="L",
                            ylabel="difference in mean winding",
                            title="Temporal mean winding minus baseline")
        for key in comparison_keys
            rows = sort([row for row in comparisons
                         if String(row.distribution) == key[1] &&
                            String(row.distribution_params) == key[2]]; by=row -> row.L)
            label = model_label(key...)
            ratios = Float64.(getproperty.(rows, :variance_ratio))
            ratio_low = Float64.(getproperty.(rows, :variance_ratio_ci_low))
            ratio_high = Float64.(getproperty.(rows, :variance_ratio_ci_high))
            plot!(ratio_plot, getproperty.(rows, :L), ratios;
                  ribbon=(ratios .- ratio_low, ratio_high .- ratios),
                  marker=:circle, label)
            differences = Float64.(getproperty.(rows, :mean_winding_difference))
            difference_low = Float64.(getproperty.(rows, :mean_winding_difference_ci_low))
            difference_high = Float64.(getproperty.(rows, :mean_winding_difference_ci_high))
            plot!(winding_plot, getproperty.(rows, :L), differences;
                  ribbon=(differences .- difference_low,
                          difference_high .- differences),
                  marker=:circle, label)
        end
        hline!(ratio_plot, [1.0]; color=:black, linestyle=:dash, label="no difference")
        hline!(winding_plot, [0.0]; color=:black, linestyle=:dash, label="no difference")
        savefig(ratio_plot, joinpath(figures_dir, "temporal_variance_ratios.png"))
        savefig(ratio_plot, joinpath(figures_dir, "temporal_variance_ratios.pdf"))
        savefig(winding_plot, joinpath(figures_dir, "temporal_mean_winding_differences.png"))
        savefig(winding_plot, joinpath(figures_dir, "temporal_mean_winding_differences.pdf"))
    end

    for name in ("summary.csv", "loglog_fits.csv", "scaling_model_comparison.csv",
                 "pointwise_ratios.csv", "local_effective_exponents.csv",
                 "temporal_baseline_comparisons.csv",
                 "temporal_direction_diagnostics.csv")
        isfile(joinpath(analysis_dir, name)) || continue
        cp(joinpath(analysis_dir, name), joinpath(output_dir, name); force=true)
    end
    pair_path = joinpath(analysis_dir, "double_dimer_pairs.csv")
    isfile(pair_path) && cp(pair_path, joinpath(output_dir, "double_dimer_pairs.csv"); force=true)
    println("Wrote Julia reports to $output_dir")
end

main()
