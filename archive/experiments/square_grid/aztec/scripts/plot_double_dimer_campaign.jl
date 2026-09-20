#!/usr/bin/env julia

using Printf

# All figures are generated as plain SVG with no plotting dependency.  The
# helper functions below share one log-scaled frame so visual encodings remain
# consistent across the difference, covariance, and decomposition plots.

function parse_arguments(arguments)
    options = Dict{String,String}(
        "analysis-dir" => joinpath(@__DIR__, "..", "output", "double_dimer_analysis"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "double_dimer_analysis"),
    )
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument in ("-h", "--help")
            println(
                "Usage: plot_double_dimer_campaign.jl " *
                "[--analysis-dir PATH] [--output-dir PATH]",
            )
            return nothing
        end
        startswith(argument, "--") || error("unexpected argument: $argument")
        key = argument[3:end]
        haskey(options, key) || error("unknown option: --$key")
        index < length(arguments) || error("missing value after --$key")
        options[key] = arguments[index + 1]
        index += 2
    end
    return (
        analysis_dir=abspath(options["analysis-dir"]),
        output_dir=abspath(options["output-dir"]),
    )
end

function read_numeric_csv(path)
    # Strict numeric-only reader for the analysis files produced in this repo.
    lines = readlines(path)
    isempty(lines) && error("empty CSV: $path")
    headers = Symbol.(split(first(lines), ','))
    rows = NamedTuple[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        values = parse.(Float64, split(line, ','))
        length(values) == length(headers) || error("malformed row in $path")
        push!(rows, NamedTuple{Tuple(headers)}(Tuple(values)))
    end
    return rows
end

function read_key_values(path)
    # variance_component_fits.txt contains only numeric key=value records.
    values = Dict{String,Float64}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        key, value = split(line, '='; limit=2)
        values[key] = parse(Float64, value)
    end
    return values
end

polyline_points(xs, ys) =
    join((@sprintf("%.3f,%.3f", x, y) for (x, y) in zip(xs, ys)), " ")

function plot_frame(io, title, description, orders, minimum_y, maximum_y)
    isempty(orders) && error("cannot plot an empty order sequence")
    maximum_y > minimum_y || error("plot y-range must be non-degenerate")
    width = 1100.0
    height = 680.0
    left, right, top, bottom = 95.0, 35.0, 70.0, 90.0
    plot_width = width - left - right
    plot_height = height - top - bottom
    minimum_log = log10(minimum(orders))
    maximum_log = log10(maximum(orders))
    x_position(order) =
        left + (log10(order) - minimum_log) / (maximum_log - minimum_log) * plot_width
    y_position(value) =
        top + (maximum_y - value) / (maximum_y - minimum_y) * plot_height

    println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $width $height\">")
    println(io, "<title>$title</title><desc>$description</desc>")
    println(io, "<rect width=\"$width\" height=\"$height\" fill=\"white\"/>")
    println(
        io,
        "<text x=\"$(width / 2)\" y=\"32\" text-anchor=\"middle\" " *
        "font-family=\"Arial,sans-serif\" font-size=\"20\" " *
        "font-weight=\"bold\">$title</text>",
    )

    for fraction in 0:5
        value = minimum_y + (maximum_y - minimum_y) * fraction / 5
        y = y_position(value)
        println(
            io,
            "<line x1=\"$left\" y1=\"$y\" x2=\"$(left + plot_width)\" " *
            "y2=\"$y\" stroke=\"#dddddd\"/>",
        )
        println(
            io,
            "<text x=\"$(left - 12)\" y=\"$(y + 5)\" text-anchor=\"end\" " *
            "font-family=\"Arial,sans-serif\" font-size=\"13\">" *
            "$(@sprintf("%.1f", value))</text>",
        )
    end
    println(
        io,
        "<line x1=\"$left\" y1=\"$top\" x2=\"$left\" " *
        "y2=\"$(top + plot_height)\" stroke=\"black\" stroke-width=\"1.5\"/>",
    )
    println(
        io,
        "<line x1=\"$left\" y1=\"$(top + plot_height)\" " *
        "x2=\"$(left + plot_width)\" y2=\"$(top + plot_height)\" " *
        "stroke=\"black\" stroke-width=\"1.5\"/>",
    )

    tick_indices = unique(
        round.(Int, range(1, length(orders), length=min(10, length(orders))))
    )
    for index in tick_indices
        order = orders[index]
        x = x_position(order)
        println(
            io,
            "<line x1=\"$x\" y1=\"$(top + plot_height)\" x2=\"$x\" " *
            "y2=\"$(top + plot_height + 6)\" stroke=\"black\"/>",
        )
        println(
            io,
            "<text x=\"$x\" y=\"$(top + plot_height + 25)\" " *
            "text-anchor=\"middle\" font-family=\"Arial,sans-serif\" " *
            "font-size=\"12\">$(Int(round(order)))</text>",
        )
    end
    println(
        io,
        "<text x=\"$(left + plot_width / 2)\" y=\"$(height - 25)\" " *
        "text-anchor=\"middle\" font-family=\"Arial,sans-serif\" " *
        "font-size=\"16\">Aztec-diamond order L (log scale)</text>",
    )
    println(
        io,
        "<text x=\"22\" y=\"$(top + plot_height / 2)\" text-anchor=\"middle\" " *
        "font-family=\"Arial,sans-serif\" font-size=\"16\" " *
        "transform=\"rotate(-90 22 $(top + plot_height / 2))\">" *
        "Height variance</text>",
    )
    return (
        x=x_position,
        y=y_position,
        left=left,
        top=top,
        plot_width=plot_width,
        plot_height=plot_height,
    )
end

function draw_series!(io, frame, orders, values; color, label, marker="circle", dash="")
    xs = frame.x.(orders)
    ys = frame.y.(values)
    dash_attribute = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
    println(
        io,
        "<polyline points=\"$(polyline_points(xs, ys))\" fill=\"none\" " *
        "stroke=\"$color\" stroke-width=\"2.5\"$dash_attribute/>",
    )
    for (x, y) in zip(xs, ys)
        if marker == "none"
            continue
        elseif marker == "square"
            println(
                io,
                "<rect x=\"$(x - 3.5)\" y=\"$(y - 3.5)\" width=\"7\" " *
                "height=\"7\" fill=\"$color\"/>",
            )
        elseif marker == "diamond"
            println(
                io,
                "<polygon points=\"$x,$(y - 5) $(x + 5),$y " *
                "$x,$(y + 5) $(x - 5),$y\" fill=\"$color\"/>",
            )
        else
            println(io, "<circle cx=\"$x\" cy=\"$y\" r=\"4\" fill=\"$color\"/>")
        end
    end
    return (color=color, label=label, dash=dash)
end

function draw_error_bars!(io, frame, orders, lows, highs; color)
    for (order, low, high) in zip(orders, lows, highs)
        x = frame.x(order)
        y_low = frame.y(low)
        y_high = frame.y(high)
        println(
            io,
            "<line x1=\"$x\" y1=\"$y_low\" x2=\"$x\" y2=\"$y_high\" " *
            "stroke=\"$color\" stroke-width=\"1.1\"/>",
        )
    end
end

function draw_legend!(io, frame, specs)
    x = frame.left + 18
    for (index, spec) in enumerate(specs)
        y = frame.top + 22 * index
        dash_attribute = isempty(spec.dash) ? "" : " stroke-dasharray=\"$(spec.dash)\""
        println(
            io,
            "<line x1=\"$x\" y1=\"$y\" x2=\"$(x + 28)\" y2=\"$y\" " *
            "stroke=\"$(spec.color)\" stroke-width=\"2.5\"$dash_attribute/>",
        )
        println(
            io,
            "<text x=\"$(x + 37)\" y=\"$(y + 5)\" " *
            "font-family=\"Arial,sans-serif\" font-size=\"13\">" *
            "$(spec.label)</text>",
        )
    end
end

function write_fit_plot(path, summary, curves)
    orders = [row.order for row in summary]
    variances = [row.variance_difference for row in summary]
    lows = [row.variance_difference_low for row in summary]
    highs = [row.variance_difference_high for row in summary]
    maximum_y = 1.08 * maximum(highs)
    open(path, "w") do io
        frame = plot_frame(
            io,
            "Double-dimer center-height difference",
            "Variance of H1-H2 with bootstrap intervals and three fitted growth curves.",
            orders,
            0.0,
            maximum_y,
        )
        draw_error_bars!(io, frame, orders, lows, highs; color="#222222")
        data_spec = draw_series!(io, frame, orders, variances; color="#222222", label="Var(H1-H2)")
        curve_orders = [row.order for row in curves]
        specs = [data_spec]
        push!(
            specs,
            draw_series!(
                io, frame, curve_orders, [row.log_fit for row in curves];
                color="#2864dc", label="a + b log L", marker="none",
            ),
        )
        push!(
            specs,
            draw_series!(
                io, frame, curve_orders, [row.log2_fit for row in curves];
                color="#e68613", label="a + b (log L)^2", marker="none",
            ),
        )
        push!(
            specs,
            draw_series!(
                io, frame, curve_orders, [row.power_fit for row in curves];
                color="#228b22", label="C (log L)^p", marker="none", dash="6,4",
            ),
        )
        draw_legend!(io, frame, specs)
        println(io, "</svg>")
    end
end

function write_covariance_fit_plot(path, summary, fits)
    # Component fits are stored separately from the primary difference fit.
    # Reconstructing curves from the reported coefficients ensures the figure
    # and machine-readable fit report cannot silently disagree.
    minimum_order = Int(round(fits["fit_min_order"]))
    fitted_rows = filter(row -> row.order >= minimum_order, summary)
    orders = [row.order for row in fitted_rows]
    covariance = [row.covariance_h1_h2 for row in fitted_rows]
    lows = [row.covariance_low for row in fitted_rows]
    highs = [row.covariance_high for row in fitted_rows]
    log_orders = log.(orders)
    prefix = "disorder_paired_covariance"
    log_curve =
        fits["$(prefix)_log_intercept"] .+ fits["$(prefix)_log_slope"] .* log_orders
    log2_curve =
        fits["$(prefix)_log2_intercept"] .+ fits["$(prefix)_log2_slope"] .* log_orders .^ 2
    power_curve =
        fits["$(prefix)_power_prefactor"] .* log_orders .^ fits["$(prefix)_power_exponent"]
    minimum_y = min(0.0, 1.05 * minimum(lows))
    maximum_y = 1.08 * maximum(vcat(highs, log_curve, log2_curve, power_curve))

    open(path, "w") do io
        frame = plot_frame(
            io,
            "Disorder-induced center-height variance",
            "Shared-environment covariance of two conditionally independent " *
            "center heights, with fitted growth curves.",
            orders,
            minimum_y,
            maximum_y,
        )
        draw_error_bars!(io, frame, orders, lows, highs; color="#222222")
        specs = NamedTuple[]
        push!(
            specs,
            draw_series!(
                io, frame, orders, covariance;
                color="#222222", label="Cov(H1,H2)", marker="diamond",
            ),
        )
        push!(specs, draw_series!(
            io, frame, orders, log_curve;
            color="#2864dc", label="a + b log L", marker="none",
        ))
        push!(specs, draw_series!(
            io, frame, orders, log2_curve;
            color="#e68613", label="a + b (log L)^2", marker="none",
        ))
        push!(specs, draw_series!(
            io, frame, orders, power_curve;
            color="#228b22", label="C (log L)^p", marker="none", dash="6,4",
        ))
        draw_legend!(io, frame, specs)
        println(io, "</svg>")
    end
end


function write_decomposition_plot(path, summary)
    orders = [row.order for row in summary]
    total = [row.total_single_variance for row in summary]
    conditional = [row.conditional_tiling_variance for row in summary]
    paired_disorder = [row.covariance_h1_h2 for row in summary]
    independent_remainder = [row.disorder_variance for row in summary]
    all_lows = vcat(
        [row.total_single_low for row in summary],
        [row.conditional_tiling_low for row in summary],
        [row.covariance_low for row in summary],
        [row.disorder_low for row in summary],
    )
    all_highs = vcat(
        [row.total_single_high for row in summary],
        [row.conditional_tiling_high for row in summary],
        [row.covariance_high for row in summary],
        [row.disorder_high for row in summary],
    )
    minimum_y = min(0.0, 1.05 * minimum(all_lows))
    maximum_y = 1.08 * maximum(all_highs)
    open(path, "w") do io
        frame = plot_frame(
            io,
            "Single-height variance decomposition",
            "Total variance split into conditional tiling and disorder-induced " *
            "components; paired covariance is the direct disorder estimate.",
            orders,
            minimum_y,
            maximum_y,
        )
        draw_error_bars!(
            io, frame, orders,
            [row.total_single_low for row in summary],
            [row.total_single_high for row in summary]; color="#2864dc",
        )
        draw_error_bars!(
            io, frame, orders,
            [row.conditional_tiling_low for row in summary],
            [row.conditional_tiling_high for row in summary]; color="#e68613",
        )
        draw_error_bars!(
            io, frame, orders,
            [row.covariance_low for row in summary],
            [row.covariance_high for row in summary]; color="#9b4dcc",
        )
        draw_error_bars!(
            io, frame, orders,
            [row.disorder_low for row in summary],
            [row.disorder_high for row in summary]; color="#228b22",
        )
        specs = NamedTuple[]
        push!(specs, draw_series!(
            io, frame, orders, total;
            color="#2864dc", label="Total Var(H)", marker="circle",
        ))
        push!(specs, draw_series!(
            io, frame, orders, conditional;
            color="#e68613", label="Tiling: Var(H1-H2)/2", marker="square",
        ))
        push!(specs, draw_series!(
            io, frame, orders, paired_disorder;
            color="#9b4dcc", label="Disorder: Cov(H1,H2)", marker="diamond",
        ))
        push!(specs, draw_series!(
            io, frame, orders, independent_remainder;
            color="#228b22", label="Independent-data remainder",
            marker="none", dash="6,4",
        ))
        draw_legend!(io, frame, specs)
        println(io, "</svg>")
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && return
    summary = read_numeric_csv(joinpath(parsed.analysis_dir, "double_dimer_summary.csv"))
    curves = read_numeric_csv(joinpath(parsed.analysis_dir, "double_dimer_fit_curves.csv"))
    component_fits =
        read_key_values(joinpath(parsed.analysis_dir, "variance_component_fits.txt"))
    mkpath(parsed.output_dir)
    fit_path = joinpath(parsed.output_dir, "double_dimer_variance_fits.svg")
    covariance_path = joinpath(parsed.output_dir, "disorder_covariance_fits.svg")
    decomposition_path = joinpath(parsed.output_dir, "variance_decomposition.svg")
    write_fit_plot(fit_path, summary, curves)
    write_covariance_fit_plot(covariance_path, summary, component_fits)
    write_decomposition_plot(decomposition_path, summary)
    println("Wrote $fit_path")
    println("Wrote $covariance_path")
    println("Wrote $decomposition_path")
end

main(ARGS)
