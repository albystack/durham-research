#!/usr/bin/env julia

using Printf

const HEADER = "order,sample_id,seed,center_row,center_column,center_height"

# Merging is also a validation boundary.  Production analyses should consume
# one canonical, sorted CSV rather than silently accepting duplicate sample IDs
# or seeds from overlapping campaigns.

function print_help()
    println("""
    Merge center-height CSV files and campaign batch directories.

    Usage:
      julia --project=aztec aztec/scripts/merge_height_batches.jl \\
        --inputs PATH1,PATH2,... --output OUTPUT.csv

    Each input may be a compact CSV or a campaign directory. Directories are
    searched recursively for batch_*.csv files. Duplicate (order, sample_id)
    keys or duplicate random seeds are rejected.
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return nothing
    options = Dict{String,String}(
        "inputs" => "",
        "output" => joinpath(
            @__DIR__,
            "..",
            "output",
            "merged_center_height_samples.csv",
        ),
    )

    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        startswith(argument, "--") || error("unexpected argument: $argument")
        key = argument[3:end]
        haskey(options, key) || error("unknown option: --$key")
        index < length(arguments) || error("missing value after --$key")
        options[key] = arguments[index + 1]
        index += 2
    end

    inputs = [
        abspath(strip(path))
        for path in split(options["inputs"], ',')
        if !isempty(strip(path))
    ]
    isempty(inputs) && error("--inputs must contain at least one file or directory")
    return (inputs=inputs, output=abspath(options["output"]))
end

function collect_csv_files(inputs)
    files = String[]
    for input in inputs
        if isfile(input)
            push!(files, input)
        elseif isdir(input)
            # Metadata and already-merged CSVs are intentionally ignored when
            # a campaign directory is supplied; only atomic batch files count.
            for (root, _, names) in walkdir(input)
                for name in names
                    startswith(name, "batch_") && endswith(name, ".csv") || continue
                    push!(files, joinpath(root, name))
                end
            end
        else
            error("input does not exist: $input")
        end
    end
    unique!(files)
    sort!(files)
    isempty(files) && error("no CSV files found")
    return files
end

function read_rows(files)
    rows = NamedTuple{
        (:order, :sample_id, :seed, :center_row, :center_column, :center_height),
        Tuple{Int,Int,UInt64,Int,Int,Int},
    }[]
    keys_seen = Set{Tuple{Int,Int}}()
    seeds_seen = Set{UInt64}()

    for path in files
        lines = readlines(path)
        isempty(lines) && error("empty CSV: $path")
        strip(first(lines)) == HEADER || error("unexpected header in $path")

        for (offset, line) in enumerate(lines[2:end])
            line_number = offset + 1
            isempty(strip(line)) && continue
            fields = split(strip(line), ',')
            length(fields) == 6 || error("malformed row $line_number in $path")
            row = (
                order=parse(Int, fields[1]),
                sample_id=parse(Int, fields[2]),
                seed=parse(UInt64, fields[3]),
                center_row=parse(Int, fields[4]),
                center_column=parse(Int, fields[5]),
                center_height=parse(Int, fields[6]),
            )

            row.order > 0 || error("non-positive order in $path line $line_number")
            row.sample_id > 0 || error("non-positive sample ID in $path line $line_number")
            row.center_row == row.order + 1 ||
                error("wrong center row in $path line $line_number")
            row.center_column == fld(row.order, 2) + 1 ||
                error("wrong center column in $path line $line_number")

            key = (row.order, row.sample_id)
            key in keys_seen && error("duplicate (order, sample_id): $key")
            row.seed in seeds_seen && error("duplicate seed: $(row.seed)")
            push!(keys_seen, key)
            push!(seeds_seen, row.seed)
            push!(rows, row)
        end
    end
    sort!(rows; by=row -> (row.order, row.sample_id))
    return rows
end

function write_rows(path, rows)
    mkpath(dirname(path))
    temporary_path = path * ".tmp"
    try
        open(temporary_path, "w") do io
            println(io, HEADER)
            for row in rows
                println(
                    io,
                    row.order, ',', row.sample_id, ',', row.seed, ',',
                    row.center_row, ',', row.center_column, ',', row.center_height,
                )
            end
        end
        # Preserve the previous merged file until the replacement is complete.
        mv(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    if isnothing(parsed)
        print_help()
        return
    end

    files = collect_csv_files(parsed.inputs)
    rows = read_rows(files)
    write_rows(parsed.output, rows)

    counts = Dict{Int,Int}()
    for row in rows
        counts[row.order] = get(counts, row.order, 0) + 1
    end
    order_summary = [
        "$(order):$(counts[order])"
        for order in sort(collect(keys(counts)))
    ]
    println("Merged $(length(rows)) independent samples from $(length(files)) files")
    println("  orders: " * join(order_summary, ", "))
    println("  output: $(parsed.output)")
end

main(ARGS)
