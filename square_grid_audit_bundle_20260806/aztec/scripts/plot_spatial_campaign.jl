#!/usr/bin/env julia

using Printf

function parse_arguments(arguments)
    options = Dict(
        "analysis-dir" => joinpath(@__DIR__, "..", "output", "spatial_analysis"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "spatial_analysis"),
    )
    index = 1
    while index <= length(arguments)
        if arguments[index] in ("-h", "--help")
            println("Usage: plot_spatial_campaign.jl [--analysis-dir PATH] [--output-dir PATH]")
            return nothing
        end
        key = replace(arguments[index], "--" => "")
        haskey(options, key) || error("unknown option: $(arguments[index])")
        index < length(arguments) || error("missing value after $(arguments[index])")
        options[key] = arguments[index + 1]
        index += 2
    end
    return (analysis_dir=abspath(options["analysis-dir"]), output_dir=abspath(options["output-dir"]))
end

function read_summary(path)
    rows = NamedTuple[]
    lines = readlines(path)
    for line in lines[2:end]
        fields = split(line, ',')
        push!(rows, (
            model=fields[1], fraction_num=parse(Int, fields[2]),
            fraction_den=parse(Int, fields[3]), order=parse(Int, fields[4]),
            separation=parse(Int, fields[5]), n=parse(Int, fields[6]),
            marginal=parse(Float64, fields[9]), marginal_low=parse(Float64, fields[10]),
            marginal_high=parse(Float64, fields[11]), conditional=parse(Float64, fields[12]),
            conditional_low=parse(Float64, fields[13]), conditional_high=parse(Float64, fields[14]),
            disorder=parse(Float64, fields[15]), disorder_low=parse(Float64, fields[16]),
            disorder_high=parse(Float64, fields[17]),
        ))
    end
    return rows
end

function read_predictions(path)
    rows = NamedTuple[]
    lines = readlines(path)
    for line in lines[2:end]
        fields = split(line, ',')
        push!(rows, (
            model=fields[1], fraction_num=parse(Int, fields[2]),
            fraction_den=parse(Int, fields[3]), component=fields[4],
            order=parse(Int, fields[5]), separation=parse(Int, fields[6]),
            observed=parse(Float64, fields[7]), log_fit=parse(Float64, fields[8]),
            quadratic_fit=parse(Float64, fields[9]),
        ))
    end
    return rows
end

function read_pooled(path)
    rows = NamedTuple[]
    lines = readlines(path)
    for line in lines[2:end]
        fields = split(line, ',')
        push!(rows, (
            model=fields[1],
            component=fields[2],
            method=fields[3],
            coefficient=parse(Float64, fields[4]),
            low=parse(Float64, fields[5]),
            high=parse(Float64, fields[6]),
            positive_fraction=parse(Float64, fields[7]),
            delta_bic=parse(Float64, fields[8]),
        ))
    end
    return rows
end

fraction_label(row) = "$(row.fraction_num)/$(row.fraction_den)"

function svg_header(io, title, description; width=1200, height=900)
    println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $width $height\">")
    println(io, "<title>$title</title><desc>$description</desc>")
    println(io, "<rect width=\"$width\" height=\"$height\" fill=\"white\"/>")
    println(io, "<text x=\"$(width/2)\" y=\"34\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"22\" font-weight=\"bold\">$title</text>")
end

function draw_panel(io, x0, y0, width, height, title, data, predictions, component)
    values = [getproperty(row, component) for row in data]
    lows = [getproperty(row, Symbol(String(component) * "_low")) for row in data]
    highs = [getproperty(row, Symbol(String(component) * "_high")) for row in data]
    all_y = vcat(values, lows, highs, [row.log_fit for row in predictions], [row.quadratic_fit for row in predictions])
    minimum_y = min(0.0, minimum(all_y))
    maximum_y = maximum(all_y)
    padding = max(0.08 * (maximum_y - minimum_y), 0.2)
    minimum_y -= padding
    maximum_y += padding
    x_values = log.([row.separation for row in data])
    min_x, max_x = extrema(x_values)
    left, right, top, bottom = 58.0, 18.0, 42.0, 48.0
    plot_width = width - left - right
    plot_height = height - top - bottom
    xp(x) = x0 + left + (x - min_x) / (max_x - min_x) * plot_width
    yp(y) = y0 + top + (maximum_y - y) / (maximum_y - minimum_y) * plot_height

    println(io, "<rect x=\"$x0\" y=\"$y0\" width=\"$width\" height=\"$height\" fill=\"#fcfcfc\" stroke=\"#cccccc\"/>")
    println(io, "<text x=\"$(x0+width/2)\" y=\"$(y0+25)\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"16\" font-weight=\"bold\">$title</text>")
    for tick in 0:4
        y = minimum_y + tick * (maximum_y - minimum_y) / 4
        py = yp(y)
        println(io, "<line x1=\"$(x0+left)\" y1=\"$py\" x2=\"$(x0+left+plot_width)\" y2=\"$py\" stroke=\"#e6e6e6\"/>")
        println(io, "<text x=\"$(x0+left-7)\" y=\"$(py+4)\" text-anchor=\"end\" font-family=\"Arial\" font-size=\"10\">$(@sprintf("%.2g", y))</text>")
    end
    println(io, "<line x1=\"$(x0+left)\" y1=\"$(y0+top)\" x2=\"$(x0+left)\" y2=\"$(y0+top+plot_height)\" stroke=\"black\"/>")
    println(io, "<line x1=\"$(x0+left)\" y1=\"$(y0+top+plot_height)\" x2=\"$(x0+left+plot_width)\" y2=\"$(y0+top+plot_height)\" stroke=\"black\"/>")
    for row in data
        x = xp(log(row.separation))
        low = yp(getproperty(row, Symbol(String(component) * "_low")))
        high = yp(getproperty(row, Symbol(String(component) * "_high")))
        y = yp(getproperty(row, component))
        println(io, "<line x1=\"$x\" y1=\"$low\" x2=\"$x\" y2=\"$high\" stroke=\"#222222\"/>")
        println(io, "<circle cx=\"$x\" cy=\"$y\" r=\"3.5\" fill=\"#222222\"/>")
    end
    for (field, color, dash) in ((:log_fit, "#2864dc", ""), (:quadratic_fit, "#e68613", "7,4"))
        points = join(("$(xp(log(row.separation))),$(yp(getproperty(row, field)))" for row in predictions), " ")
        dash_attr = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
        println(io, "<polyline points=\"$points\" fill=\"none\" stroke=\"$color\" stroke-width=\"2.2\"$dash_attr/>")
    end
    println(io, "<text x=\"$(x0+left+plot_width/2)\" y=\"$(y0+height-10)\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"12\">log separation</text>")
end

function write_four_panel(path, summary, predictions; model, component, title)
    fractions = sort(unique((row.fraction_num, row.fraction_den) for row in summary if row.model == model); by=f -> f[1]/f[2])
    length(fractions) == 4 || error("four-panel plot expects four fractions")
    open(path, "w") do io
        svg_header(io, title, "Bootstrap intervals with nested log and log-plus-log-squared fits.")
        for (index, fraction) in enumerate(fractions)
            column = (index - 1) % 2
            row_index = (index - 1) ÷ 2
            x0 = 35 + column * 575
            y0 = 65 + row_index * 400
            data = sort([
                row for row in summary
                if row.model == model && (row.fraction_num, row.fraction_den) == fraction
            ]; by=row -> row.order)
            curves = sort([
                row for row in predictions
                if row.model == model && row.component == String(component) &&
                   (row.fraction_num, row.fraction_den) == fraction
            ]; by=row -> row.order)
            draw_panel(io, x0, y0, 540, 365, "separation = $(fraction[1])/$(fraction[2]) of L", data, curves, component)
        end
        println(io, "<line x1=\"420\" y1=\"875\" x2=\"450\" y2=\"875\" stroke=\"#2864dc\" stroke-width=\"2.2\"/><text x=\"458\" y=\"879\" font-family=\"Arial\" font-size=\"12\">a + b log r</text>")
        println(io, "<line x1=\"590\" y1=\"875\" x2=\"620\" y2=\"875\" stroke=\"#e68613\" stroke-width=\"2.2\" stroke-dasharray=\"7,4\"/><text x=\"628\" y=\"879\" font-family=\"Arial\" font-size=\"12\">a + b log r + c(log r)²</text>")
        println(io, "</svg>")
    end
end

function write_decomposition(path, summary)
    fraction = (1, 4)
    gamma = sort([row for row in summary if row.model == "gamma" && (row.fraction_num, row.fraction_den) == fraction]; by=row -> row.order)
    uniform = sort([row for row in summary if row.model == "uniform" && (row.fraction_num, row.fraction_den) == fraction]; by=row -> row.order)
    all_values = vcat(
        [row.marginal_high for row in gamma], [row.conditional_high for row in gamma],
        [row.disorder_high for row in gamma], [row.marginal_high for row in uniform],
    )
    max_y = 1.08 * maximum(all_values)
    min_x, max_x = extrema(log.([row.separation for row in gamma]))
    xp(r) = 100 + (log(r) - min_x) / (max_x - min_x) * 1030
    yp(y) = 70 + (max_y - y) / max_y * 520
    open(path, "w") do io
        svg_header(io, "Spatial variance decomposition at separation L/4", "Gamma thermal and disorder terms compared with a uniform control."; width=1200, height=700)
        for tick in 0:5
            y = tick * max_y / 5
            println(io, "<line x1=\"100\" y1=\"$(yp(y))\" x2=\"1130\" y2=\"$(yp(y))\" stroke=\"#e2e2e2\"/>")
            println(io, "<text x=\"88\" y=\"$(yp(y)+4)\" text-anchor=\"end\" font-family=\"Arial\" font-size=\"12\">$(@sprintf("%.2g", y))</text>")
        end
        println(io, "<line x1=\"100\" y1=\"70\" x2=\"100\" y2=\"590\" stroke=\"black\"/><line x1=\"100\" y1=\"590\" x2=\"1130\" y2=\"590\" stroke=\"black\"/>")
        series = (
            (rows=gamma, field=:marginal, color="#222222", label="Gamma total"),
            (rows=gamma, field=:conditional, color="#2864dc", label="Gamma conditional"),
            (rows=gamma, field=:disorder, color="#e68613", label="Gamma disorder covariance"),
            (rows=uniform, field=:marginal, color="#228b22", label="Uniform control"),
        )
        for item in series
            points = join(("$(xp(row.separation)),$(yp(getproperty(row, item.field)))" for row in item.rows), " ")
            println(io, "<polyline points=\"$points\" fill=\"none\" stroke=\"$(item.color)\" stroke-width=\"2.5\"/>")
            for row in item.rows
                println(io, "<circle cx=\"$(xp(row.separation))\" cy=\"$(yp(getproperty(row, item.field)))\" r=\"3.5\" fill=\"$(item.color)\"/>")
            end
        end
        println(io, "<text x=\"615\" y=\"630\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"15\">separation r (log scale)</text>")
        for (index, item) in enumerate(series)
            x = 160 + 250 * ((index - 1) % 4)
            println(io, "<line x1=\"$x\" y1=\"665\" x2=\"$(x+28)\" y2=\"665\" stroke=\"$(item.color)\" stroke-width=\"2.5\"/><text x=\"$(x+35)\" y=\"669\" font-family=\"Arial\" font-size=\"12\">$(item.label)</text>")
        end
        println(io, "</svg>")
    end
end

function write_pooled_coefficients(path, rows)
    selected = [
        row for row in rows
        if (row.model == "gamma" && row.component in ("disorder", "conditional")) ||
           (row.model == "uniform" && row.component == "marginal")
    ]
    minimum_x = min(-0.25, minimum(row.low for row in selected))
    maximum_x = max(0.25, maximum(row.high for row in selected))
    padding = 0.12 * (maximum_x - minimum_x)
    minimum_x -= padding
    maximum_x += padding
    x_position(value) = 250 + (value - minimum_x) / (maximum_x - minimum_x) * 870
    labels = (
        (model="gamma", component="disorder", text="Gamma disorder covariance"),
        (model="gamma", component="conditional", text="Gamma conditional component"),
        (model="uniform", component="marginal", text="Uniform marginal control"),
    )
    open(path, "w") do io
        svg_header(
            io,
            "Common log-squared coefficient across four separations",
            "Environment-clustered pooled coefficient with 95 percent bootstrap intervals.";
            width=1200,
            height=520,
        )
        zero_x = x_position(0)
        println(io, "<line x1=\"$zero_x\" y1=\"85\" x2=\"$zero_x\" y2=\"410\" stroke=\"#555555\" stroke-width=\"1.5\" stroke-dasharray=\"5,4\"/>")
        for tick in 0:5
            value = minimum_x + tick * (maximum_x - minimum_x) / 5
            x = x_position(value)
            println(io, "<line x1=\"$x\" y1=\"410\" x2=\"$x\" y2=\"417\" stroke=\"black\"/>")
            println(io, "<text x=\"$x\" y=\"437\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"12\">$(@sprintf("%.2g", value))</text>")
        end
        println(io, "<line x1=\"250\" y1=\"410\" x2=\"1120\" y2=\"410\" stroke=\"black\"/>")
        for (label_index, label) in enumerate(labels)
            base_y = 125 + (label_index - 1) * 125
            println(io, "<text x=\"230\" y=\"$(base_y+5)\" text-anchor=\"end\" font-family=\"Arial\" font-size=\"15\">$(label.text)</text>")
            for (method_index, method) in enumerate(("unweighted", "weighted"))
                row = only([
                    value for value in selected
                    if value.model == label.model && value.component == label.component &&
                       value.method == method
                ])
                y = base_y + (method_index == 1 ? -15 : 15)
                color = method == "unweighted" ? "#2864dc" : "#e68613"
                println(io, "<line x1=\"$(x_position(row.low))\" y1=\"$y\" x2=\"$(x_position(row.high))\" y2=\"$y\" stroke=\"$color\" stroke-width=\"3\"/>")
                println(io, "<line x1=\"$(x_position(row.low))\" y1=\"$(y-6)\" x2=\"$(x_position(row.low))\" y2=\"$(y+6)\" stroke=\"$color\"/>")
                println(io, "<line x1=\"$(x_position(row.high))\" y1=\"$(y-6)\" x2=\"$(x_position(row.high))\" y2=\"$(y+6)\" stroke=\"$color\"/>")
                println(io, "<circle cx=\"$(x_position(row.coefficient))\" cy=\"$y\" r=\"6\" fill=\"$color\"/>")
            end
        end
        println(io, "<text x=\"685\" y=\"475\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"15\">coefficient c in a + b log(r) + c(log(r))²</text>")
        println(io, "<circle cx=\"470\" cy=\"502\" r=\"5\" fill=\"#2864dc\"/><text x=\"482\" y=\"507\" font-family=\"Arial\" font-size=\"12\">unweighted</text>")
        println(io, "<circle cx=\"610\" cy=\"502\" r=\"5\" fill=\"#e68613\"/><text x=\"622\" y=\"507\" font-family=\"Arial\" font-size=\"12\">inverse-variance weighted</text>")
        println(io, "</svg>")
    end
end

function maybe_write_png(svg_path)
    # macOS ships `sips`, which can rasterise the dependency-free SVG figures.
    # CI and other platforms still retain the publication-quality SVG files.
    sips = Sys.which("sips")
    isnothing(sips) && return nothing
    png_path = replace(svg_path, r"\.svg$" => ".png")
    command = `$sips -s format png $svg_path --out $png_path`
    success(pipeline(command; stdout=devnull, stderr=devnull)) ||
        @warn "could not rasterise SVG" svg_path
    return png_path
end

function main(arguments)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && return
    summary = read_summary(joinpath(parsed.analysis_dir, "spatial_summary.csv"))
    predictions = read_predictions(joinpath(parsed.analysis_dir, "spatial_fit_curves.csv"))
    pooled = read_pooled(joinpath(parsed.analysis_dir, "spatial_pooled_model_comparison.csv"))
    mkpath(parsed.output_dir)
    gamma_path = joinpath(parsed.output_dir, "gamma_disorder_spatial_fits.svg")
    uniform_path = joinpath(parsed.output_dir, "uniform_control_spatial_fits.svg")
    decomposition_path = joinpath(parsed.output_dir, "spatial_variance_decomposition.svg")
    pooled_path = joinpath(parsed.output_dir, "pooled_log2_coefficients.svg")
    write_four_panel(
        gamma_path,
        summary,
        predictions;
        model="gamma",
        component=:disorder,
        title="Gamma disorder-induced spatial height fluctuations",
    )
    write_four_panel(
        uniform_path,
        summary,
        predictions;
        model="uniform",
        component=:marginal,
        title="Uniform Aztec-diamond spatial height control",
    )
    write_decomposition(decomposition_path, summary)
    write_pooled_coefficients(pooled_path, pooled)
    maybe_write_png(gamma_path)
    maybe_write_png(uniform_path)
    maybe_write_png(decomposition_path)
    maybe_write_png(pooled_path)
    println("Spatial figures written to $(parsed.output_dir)")
end

main(ARGS)
