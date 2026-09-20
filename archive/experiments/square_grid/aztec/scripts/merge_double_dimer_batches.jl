#!/usr/bin/env julia

const HEADER =
    "order,sample_id,seed,center_row,center_column,height_1,height_2,height_difference"

# Treat merge time as a strict data-integrity check.  In particular, the
# stored difference is recomputed from the two raw heights for every row.

function parse_arguments(arguments)
    options = Dict{String,String}(
        "inputs" => "",
        "output" => joinpath(@__DIR__, "..", "output", "double_dimer_pairs.csv"),
    )
    if any(argument -> argument in ("-h", "--help"), arguments)
        println("""
        Usage:
          julia --project=aztec aztec/scripts/merge_double_dimer_batches.jl \\
            --inputs PATH1,PATH2,... --output OUTPUT.csv
        """)
        return nothing
    end
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
    isempty(inputs) && error("--inputs must contain at least one path")
    return (inputs=inputs, output=abspath(options["output"]))
end

function collect_files(inputs)
    files = String[]
    for input in inputs
        if isfile(input)
            push!(files, input)
        elseif isdir(input)
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
    isempty(files) && error("no paired CSV files found")
    return files
end


function load_rows(files)
    rows = NamedTuple{
        (:order, :sample_id, :seed, :center_row, :center_column,
         :height_1, :height_2, :height_difference),
        Tuple{Int,Int,UInt64,Int,Int,Int,Int,Int},
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
            length(fields) == 8 || error("malformed row $line_number in $path")
            row = (
                order=parse(Int, fields[1]),
                sample_id=parse(Int, fields[2]),
                seed=parse(UInt64, fields[3]),
                center_row=parse(Int, fields[4]),
                center_column=parse(Int, fields[5]),
                height_1=parse(Int, fields[6]),
                height_2=parse(Int, fields[7]),
                height_difference=parse(Int, fields[8]),
            )
            row.order > 0 || error("non-positive order in $path line $line_number")
            row.sample_id > 0 ||
                error("non-positive sample ID in $path line $line_number")
            row.center_row == row.order + 1 ||
                error("wrong center row in $path line $line_number")
            row.center_column == fld(row.order, 2) + 1 ||
                error("wrong center column in $path line $line_number")
            row.height_difference == row.height_1 - row.height_2 ||
                error("wrong height difference in $path line $line_number")
            key = (row.order, row.sample_id)
            key in keys_seen && error("duplicate paired key $key")
            row.seed in seeds_seen && error("duplicate paired seed $(row.seed)")
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
                    row.center_row, ',', row.center_column, ',', row.height_1, ',',
                    row.height_2, ',', row.height_difference,
                )
            end
        end
        # Atomic replacement makes a failed merge recoverable.
        mv(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && return
    files = collect_files(parsed.inputs)
    rows = load_rows(files)
    write_rows(parsed.output, rows)
    counts = Dict{Int,Int}()
    for row in rows
        counts[row.order] = get(counts, row.order, 0) + 1
    end
    summary = ["$(order):$(counts[order])" for order in sort(collect(keys(counts)))]
    println("Merged $(length(rows)) double-dimer pairs from $(length(files)) files")
    println("  orders: " * join(summary, ", "))
    println("  output: $(parsed.output)")
end

main(ARGS)
