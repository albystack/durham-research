#!/usr/bin/env julia

using Printf

const HEADER =
    "model,order,sample_id,seed,fraction_num,fraction_den,separation," *
    "left_column,right_column,increment_1,increment_2,difference"

function print_help()
    println("""
    Merge validated square-grid batch CSVs.

    Usage:
      julia --project=aztec aztec/scripts/merge_square_grid_batches.jl \\
        --inputs DIR[,DIR...] --output PATH
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return nothing
    options = Dict("inputs" => "", "output" => "square_grid_pairs.csv")
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
    isempty(options["inputs"]) && error("--inputs is required")
    return (
        inputs=[abspath(strip(path)) for path in split(options["inputs"], ',')],
        output=abspath(options["output"]),
    )
end

function collect_files(paths)
    files = String[]
    for path in paths
        if isfile(path)
            push!(files, path)
        elseif isdir(path)
            for (root, _, names) in walkdir(path), name in names
                startswith(name, "batch_") && endswith(name, ".csv") || continue
                push!(files, joinpath(root, name))
            end
        else
            error("input path does not exist: $path")
        end
    end
    unique!(files)
    sort!(files)
    isempty(files) && error("no batch CSVs found")
    return files
end

function load_rows(files)
    rows = NamedTuple[]
    seen = Set{Tuple{String,Int,Int,Int,Int}}()
    seed_by_sample = Dict{Tuple{String,Int,Int},UInt64}()
    for path in files
        lines = readlines(path)
        isempty(lines) && error("empty batch: $path")
        strip(first(lines)) == HEADER || error("unexpected header: $path")
        for (offset, line) in enumerate(lines[2:end])
            isempty(strip(line)) && continue
            fields = split(strip(line), ',')
            length(fields) == 12 || error("malformed row $(offset + 1) in $path")
            row = (
                model=fields[1],
                order=parse(Int, fields[2]),
                sample_id=parse(Int, fields[3]),
                seed=parse(UInt64, fields[4]),
                fraction_num=parse(Int, fields[5]),
                fraction_den=parse(Int, fields[6]),
                separation=parse(Int, fields[7]),
                left_column=parse(Int, fields[8]),
                right_column=parse(Int, fields[9]),
                increment_1=parse(Int, fields[10]),
                increment_2=parse(Int, fields[11]),
                difference=parse(Int, fields[12]),
            )
            row.difference == row.increment_1 - row.increment_2 ||
                error("incorrect difference in $path")
            row.right_column - row.left_column == row.separation ||
                error("incorrect separation in $path")
            key = (
                row.model, row.order, row.sample_id,
                row.fraction_num, row.fraction_den)
            key in seen && error("duplicate square-grid observation $key")
            push!(seen, key)
            sample_key = (row.model, row.order, row.sample_id)
            if haskey(seed_by_sample, sample_key)
                seed_by_sample[sample_key] == row.seed ||
                    error("inconsistent seed for $sample_key")
            else
                seed_by_sample[sample_key] = row.seed
            end
            push!(rows, row)
        end
    end
    sort!(rows; by=row -> (
        row.model, row.order, row.sample_id,
        row.fraction_num / row.fraction_den))
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
                    row.model, ',', row.order, ',', row.sample_id, ',', row.seed, ',',
                    row.fraction_num, ',', row.fraction_den, ',', row.separation, ',',
                    row.left_column, ',', row.right_column, ',', row.increment_1, ',',
                    row.increment_2, ',', row.difference,
                )
            end
        end
        mv(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && (print_help(); return)
    files = collect_files(parsed.inputs)
    rows = load_rows(files)
    write_rows(parsed.output, rows)
    println("Merged $(length(rows)) square-grid spatial rows from $(length(files)) files")
    println("  output: $(parsed.output)")
end

main(ARGS)
