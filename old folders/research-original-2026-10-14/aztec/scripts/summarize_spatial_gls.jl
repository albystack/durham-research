#!/usr/bin/env julia

using Printf

function print_help()
    println("""
    Consolidate block-GLS outputs from several spatial analyses.

    Usage:
      julia summarize_spatial_gls.jl \
        --entry LABEL=ANALYSIS_DIR [--entry ...] --output-dir PATH
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && (print_help(); return nothing)
    entries = Pair{String,String}[]
    output_dir = ""
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        index < length(arguments) || error("missing value after $argument")
        value = arguments[index + 1]
        if argument == "--entry"
            fields = split(value, '='; limit=2)
            length(fields) == 2 || error("--entry must have LABEL=ANALYSIS_DIR form")
            isempty(fields[1]) && error("entry label cannot be empty")
            occursin(',', fields[1]) && error("entry label cannot contain a comma")
            push!(entries, fields[1] => abspath(fields[2]))
        elseif argument == "--output-dir"
            output_dir = abspath(value)
        else
            error("unknown option: $argument")
        end
        index += 2
    end
    isempty(entries) && error("at least one --entry is required")
    length(unique(first.(entries))) == length(entries) || error("entry labels must be unique")
    isempty(output_dir) && error("--output-dir is required")
    return (entries=entries, output_dir=output_dir)
end

function read_table(label, directory)
    path = joinpath(directory, "spatial_pooled_gls_comparison.csv")
    isfile(path) || error("missing block-GLS table: $path")
    lines = readlines(path)
    length(lines) >= 2 || error("empty block-GLS table: $path")
    header = split(first(lines), ',')
    expected = [
        "model", "component", "method", "common_log2_coefficient",
        "coefficient_low", "coefficient_high", "coefficient_positive_fraction",
        "delta_bic_log_minus_quadratic", "delta_bic_low", "delta_bic_high",
        "quadratic_bic_win_fraction", "log_holdout_rmse", "quadratic_holdout_rmse",
    ]
    header == expected || error("unexpected block-GLS schema in $path")
    return [
        begin
            fields = split(line, ',')
            length(fields) == length(header) || error("malformed row in $path")
            NamedTuple{Tuple(Symbol.(header))}(Tuple(fields))
        end
        for line in lines[2:end] if !isempty(strip(line))
    ], label
end

function main(arguments)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && return
    collected = NamedTuple[]
    for (label, directory) in parsed.entries
        rows, _ = read_table(label, directory)
        append!(collected, [(campaign=label, row...) for row in rows])
    end
    mkpath(parsed.output_dir)
    csv_path = joinpath(parsed.output_dir, "spatial_gls_campaign_comparison.csv")
    open(csv_path, "w") do io
        println(io, "campaign," * join(string.(keys(first(collected))[2:end]), ','))
        for row in collected
            println(io, join(values(row), ','))
        end
    end
    report_path = joinpath(parsed.output_dir, "spatial_gls_campaign_comparison.md")
    open(report_path, "w") do io
        println(io, "# Spatial block-GLS campaign comparison\n")
        println(io, "The table below selects the disorder covariance from each analysis. ")
        println(io, "Intervals and model comparisons come from joint environment bootstraps.\n")
        println(io, "| Campaign | min-order tag | c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / log2 |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|")
        for row in collected
            row.model == "gamma" && row.component == "disorder" || continue
            cutoff_match = match(r"min(\d+)", row.campaign)
            cutoff = isnothing(cutoff_match) ? "?" : cutoff_match.captures[1]
            @printf(
                io,
                "| %s | %s | %.4g | [%.4g, %.4g] | %.4f | %.3f | %.3f / %.3f |\n",
                row.campaign,
                cutoff,
                parse(Float64, row.common_log2_coefficient),
                parse(Float64, row.coefficient_low),
                parse(Float64, row.coefficient_high),
                parse(Float64, row.coefficient_positive_fraction),
                parse(Float64, row.delta_bic_log_minus_quadratic),
                parse(Float64, row.log_holdout_rmse),
                parse(Float64, row.quadratic_holdout_rmse),
            )
        end
        println(io, "\nThese are finite-size numerical comparisons, not asymptotic proofs.")
    end
    println("Wrote $csv_path")
    println("Wrote $report_path")
end

main(ARGS)
