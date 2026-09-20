#!/usr/bin/env julia

using Printf

# Dependency-free plotting keeps the retained analysis reproducible with a
# stock Julia installation.  The output is editable vector SVG; PNG files in
# results/ are rendered derivatives for convenient previewing.

function print_help()
    println("""
    Plot the single-height summary and fitted curves as SVG.

    Usage:
      julia --project=aztec aztec/scripts/plot_height_campaign.jl [options]

    Options:
      --analysis-dir PATH  directory containing height_summary.csv and curves
      --output-dir PATH    destination directory
      -h, --help           show this message
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return nothing
    options = Dict{String,String}(
        "analysis-dir" => joinpath(@__DIR__, "..", "output", "gamma_height_analysis"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "gamma_height_analysis"),
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
    return (
        analysis_dir=abspath(options["analysis-dir"]),
        output_dir=abspath(options["output-dir"]),
    )
end

function read_numeric_csv(path)
    # Analysis CSVs are deliberately numeric-only, so a tiny strict reader is
    # sufficient and avoids adding a dataframe/CSV dependency just for plots.
    lines = readlines(path)
    isempty(lines) && error("empty CSV: $path")
    headers = Symbol.(split(first(lines), ','))
    rows = NamedTuple[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        values = parse.(Float64, split(line, ','))
        length(values) == length(headers) || error("malformed CSV row in $path")
        push!(rows, NamedTuple{Tuple(headers)}(Tuple(values)))
    end
    return rows
end

xml_escape(text) = replace(
    string(text),
    '&' => "&amp;",
    '<' => "&lt;",
    '>' => "&gt;",
    '"' => "&quot;",
)

function polyline_points(xs, ys)
    return join((@sprintf("%.3f,%.3f", x, y) for (x, y) in zip(xs, ys)), " ")
end

function write_plot(path, summary, curves)
    width = 1100.0
    height = 680.0
    left = 95.0
    right = 35.0
    top = 70.0
    bottom = 85.0
    plot_width = width - left - right
    plot_height = height - top - bottom

    isempty(summary) && error("height summary is empty")
    isempty(curves) && error("height fit curve is empty")
    orders = [row.order for row in summary]
    variances = [row.variance_height for row in summary]
    lows = [row.variance_bootstrap_low for row in summary]
    highs = [row.variance_bootstrap_high for row in summary]
    minimum_log_order = log10(minimum(orders))
    maximum_log_order = log10(maximum(orders))
    maximum_y = 1.08 * maximum(highs)
    maximum_y > 0 || error("variance plot requires positive values")

    x_position(order) =
        left +
        (log10(order) - minimum_log_order) /
        (maximum_log_order - minimum_log_order) * plot_width
    y_position(value) = top + (1 - value / maximum_y) * plot_height

    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        println(
            io,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" ",
            "viewBox=\"0 0 $width $height\">",
        )
        println(
            io,
            "<title>Gamma-disordered Aztec diamond: center-height fluctuations</title>",
            "<desc>Monte Carlo center-height variances with bootstrap intervals and ",
            "ordinary-log, squared-log, and free-power fitted curves through order 1300.</desc>",
        )
        println(io, "<rect width=\"$width\" height=\"$height\" fill=\"white\"/>")
        println(
            io,
            "<text x=\"$(width / 2)\" y=\"32\" text-anchor=\"middle\" ",
            "font-family=\"Arial,sans-serif\" font-size=\"20\" font-weight=\"bold\">",
            "Gamma-disordered Aztec diamond: center-height fluctuations</text>",
        )

        for fraction in 0:5
            value = maximum_y * fraction / 5
            y = y_position(value)
            println(
                io,
                "<line x1=\"$left\" y1=\"$y\" x2=\"$(left + plot_width)\" y2=\"$y\" ",
                "stroke=\"#dddddd\" stroke-width=\"1\"/>",
            )
            println(
                io,
                "<text x=\"$(left - 12)\" y=\"$(y + 5)\" text-anchor=\"end\" ",
                "font-family=\"Arial,sans-serif\" font-size=\"13\">",
                @sprintf("%.1f", value),
                "</text>",
            )
        end

        println(
            io,
            "<line x1=\"$left\" y1=\"$top\" x2=\"$left\" y2=\"$(top + plot_height)\" ",
            "stroke=\"black\" stroke-width=\"1.5\"/>",
        )
        println(
            io,
            "<line x1=\"$left\" y1=\"$(top + plot_height)\" ",
            "x2=\"$(left + plot_width)\" y2=\"$(top + plot_height)\" ",
            "stroke=\"black\" stroke-width=\"1.5\"/>",
        )

        # Label at most ten logarithmically positioned observed orders; all
        # observations are still plotted even when their tick label is omitted.
        tick_indices = unique(
            round.(Int, range(1, length(orders), length=min(10, length(orders))))
        )
        for index in tick_indices
            order = orders[index]
            x = x_position(order)
            println(
                io,
                "<line x1=\"$x\" y1=\"$(top + plot_height)\" x2=\"$x\" ",
                "y2=\"$(top + plot_height + 6)\" stroke=\"black\"/>",
            )
            println(
                io,
                "<text x=\"$x\" y=\"$(top + plot_height + 25)\" text-anchor=\"middle\" ",
                "font-family=\"Arial,sans-serif\" font-size=\"12\">",
                Int(round(order)),
                "</text>",
            )
        end

        curve_specs = (
            (:log_fit, "#2864dc", "a + b log L", ""),
            (:log2_fit, "#e68613", "a + b (log L)^2", ""),
            (:power_fit, "#228b22", "C (log L)^p", "6,4"),
        )
        curve_orders = [row.order for row in curves]
        for (field, color, _, dash) in curve_specs
            curve_values = [getproperty(row, field) for row in curves]
            xs = x_position.(curve_orders)
            ys = y_position.(curve_values)
            dash_attribute = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
            println(
                io,
                "<polyline points=\"$(polyline_points(xs, ys))\" fill=\"none\" ",
                "stroke=\"$color\" stroke-width=\"2.5\"$dash_attribute/>",
            )
        end

        for (order, variance, low, high) in zip(orders, variances, lows, highs)
            x = x_position(order)
            y = y_position(variance)
            y_low = y_position(low)
            y_high = y_position(high)
            println(
                io,
                "<line x1=\"$x\" y1=\"$y_low\" x2=\"$x\" y2=\"$y_high\" ",
                "stroke=\"black\" stroke-width=\"1.2\"/>",
            )
            println(io, "<line x1=\"$(x - 4)\" y1=\"$y_low\" x2=\"$(x + 4)\" ",
                    "y2=\"$y_low\" stroke=\"black\"/>")
            println(io, "<line x1=\"$(x - 4)\" y1=\"$y_high\" x2=\"$(x + 4)\" ",
                    "y2=\"$y_high\" stroke=\"black\"/>")
            println(io, "<circle cx=\"$x\" cy=\"$y\" r=\"4\" fill=\"black\"/>")
        end

        println(
            io,
            "<text x=\"$(left + plot_width / 2)\" y=\"$(height - 25)\" ",
            "text-anchor=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"16\">",
            "Aztec-diamond order L (log scale)</text>",
        )
        println(
            io,
            "<text x=\"22\" y=\"$(top + plot_height / 2)\" text-anchor=\"middle\" ",
            "font-family=\"Arial,sans-serif\" font-size=\"16\" ",
            "transform=\"rotate(-90 22 $(top + plot_height / 2))\">",
            "Variance of center height</text>",
        )

        legend_x = left + 18
        legend_y = top + 18
        for (index, (_, color, label, dash)) in enumerate(curve_specs)
            y = legend_y + 24 * index
            dash_attribute = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
            println(
                io,
                "<line x1=\"$legend_x\" y1=\"$y\" x2=\"$(legend_x + 28)\" y2=\"$y\" ",
                "stroke=\"$color\" stroke-width=\"2.5\"$dash_attribute/>",
            )
            println(
                io,
                "<text x=\"$(legend_x + 37)\" y=\"$(y + 5)\" ",
                "font-family=\"Arial,sans-serif\" font-size=\"13\">",
                xml_escape(label),
                "</text>",
            )
        end
        point_y = legend_y
        println(io, "<circle cx=\"$(legend_x + 14)\" cy=\"$point_y\" r=\"4\" fill=\"black\"/>")
        println(
            io,
            "<text x=\"$(legend_x + 37)\" y=\"$(point_y + 5)\" ",
            "font-family=\"Arial,sans-serif\" font-size=\"13\">Monte Carlo variance</text>",
        )
        println(io, "</svg>")
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    if isnothing(parsed)
        print_help()
        return
    end
    summary = read_numeric_csv(joinpath(parsed.analysis_dir, "height_summary.csv"))
    curves = read_numeric_csv(joinpath(parsed.analysis_dir, "height_fit_curves.csv"))
    mkpath(parsed.output_dir)
    output_path = joinpath(parsed.output_dir, "center_height_variance.svg")
    write_plot(output_path, summary, curves)
    println("Wrote $output_path")
end

main(ARGS)
