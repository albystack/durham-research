#!/usr/bin/env julia

using LERWResearch

function main(args=ARGS)::Int
    options = LERWResearch.parse_cli(args)
    roots = split(LERWResearch.require_option(options, "results-dirs"), ',')
    expected_rows = parse.(Int, split(LERWResearch.require_option(options, "expected-rows"), ','))
    length(roots) == length(expected_rows) ||
        throw(ArgumentError("--results-dirs and --expected-rows must have equal lengths"))

    allow_incomplete = haskey(options, "allow-incomplete")
    raw = NamedTuple[]
    validation = NamedTuple[]

    for (root, expected) in zip(roots, expected_rows)
        isdir(root) || throw(ArgumentError("results directory does not exist: $root"))
        files = String[]
        for (directory, _, names) in walkdir(root), name in names
            startswith(name, "batch_") && endswith(name, ".csv") &&
                push!(files, joinpath(directory, name))
        end
        for path in sort!(files)
            rows = LERWResearch.csv_rows(path)
            ok_rows = [row for row in rows if row.status == "ok"]
            complete = length(rows) == expected && length(ok_rows) == expected
            push!(validation, (
                results_root=root, path=path, expected_rows=expected,
                row_count=length(rows), ok_count=length(ok_rows), complete=complete,
            ))
            (complete || allow_incomplete) && append!(raw, ok_rows)
        end
    end

    incomplete = [row.path for row in validation if !row.complete]
    !allow_incomplete && !isempty(incomplete) &&
        throw(ErrorException("incomplete result files: $(join(incomplete[1:min(end, 10)], ", "))"))
    isempty(raw) && throw(ErrorException("no valid observations found"))
    LERWResearch.assert_no_duplicates(raw)

    output_dir = get(options, "output-dir", "analysis_combined")
    bootstrap_reps = parse(Int, get(options, "bootstrap-reps", "0"))
    bootstrap_seed = parse(Int, get(options, "bootstrap-seed", "20260623"))
    fit_min_L = haskey(options, "fit-min-L") ? parse(Int, options["fit-min-L"]) : nothing

    summary = LERWResearch.make_summary(raw)
    annealed = LERWResearch.fit_exponents(summary, :annealed; min_L=fit_min_L)
    quenched = LERWResearch.fit_exponents(summary, :quenched; min_L=fit_min_L)
    LERWResearch.bootstrap_intervals!(annealed, raw, summary, :annealed,
                                      bootstrap_reps, UInt64(bootstrap_seed))
    LERWResearch.bootstrap_intervals!(quenched, raw, summary, :quenched,
                                      bootstrap_reps, UInt64(bootstrap_seed + 7919))

    LERWResearch.write_table(joinpath(output_dir, "validation.csv"), validation)
    LERWResearch.write_table(joinpath(output_dir, "combined_raw.csv"), raw)
    LERWResearch.write_table(joinpath(output_dir, "summary.csv"), summary)
    LERWResearch.write_table(joinpath(output_dir, "loglog_fits.csv"), vcat(annealed, quenched))
    LERWResearch.write_table(joinpath(output_dir, "scaling_model_comparison.csv"),
                             LERWResearch.fit_scaling_models(summary))
    LERWResearch.write_table(joinpath(output_dir, "pointwise_ratios.csv"),
                             LERWResearch.pointwise_rows(summary))
    LERWResearch.write_table(joinpath(output_dir, "local_effective_exponents.csv"),
                             LERWResearch.local_exponent_rows(summary))

    println("Analysed $(length(raw)) observations from $(length(validation)) batch files")
    println("Output: $output_dir")
    return 0
end

exit(main())
