const MODEL_NAMES = Set([
    "baseline", "symmetric", "gamma", "gamma_edges", "exponential", "exp",
    "exp_edges", "lognormal", "lognormal_edges", "pareto", "pareto_edges",
    "uniform", "uniform_edges", "beta", "beta_edges", "weibull",
    "weibull_edges", "inverse_gamma", "inverse_gamma_edges", "bernoulli",
    "bernoulli_edges", "triangular", "triangular_edges",
])

"Canonical compact JSON used in metadata and seeds."
function canonical_json(params)::String
    isempty(params) && return "{}"
    keys_sorted = sort!(String.(collect(keys(params))))
    pieces = String[]
    for key in keys_sorted
        value = haskey(params, key) ? params[key] : params[Symbol(key)]
        push!(pieces, string(JSON3.write(key), ":", JSON3.write(value)))
    end
    return string("{", join(pieces, ","), "}")
end

"A platform-independent 64-bit seed derived from SHA-256."
function stable_seed(parts...)::UInt64
    digest = sha256(join(string.(parts), ":"))
    value = UInt64(0)
    @inbounds for byte in digest[1:8]
        value = (value << 8) | UInt64(byte)
    end
    return value
end

function parse_params(text)::Dict{String,Float64}
    (text === nothing || isempty(strip(String(text)))) && return Dict{String,Float64}()
    parsed = JSON3.read(String(text))
    parsed isa JSON3.Object || throw(ArgumentError("distribution_params must be a JSON object"))
    return Dict(String(key) => Float64(value) for (key, value) in pairs(parsed))
end

function require_parameter(params::AbstractDict, aliases::Tuple, canonical::String)::Float64
    for name in aliases
        haskey(params, name) && return Float64(params[name])
    end
    throw(ArgumentError("distribution needs parameter `$canonical`"))
end

"Return `(canonical_model, parameter)` and validate its parameter range."
function distribution_spec(distribution::AbstractString, params::AbstractDict)
    name = strip(String(distribution))
    name in MODEL_NAMES || throw(ArgumentError("unknown distribution: $name"))

    model, parameter = if name in ("baseline", "symmetric")
        isempty(params) || throw(ArgumentError("$name does not take parameters"))
        (:baseline, nothing)
    elseif name in ("exponential", "exp", "exp_edges")
        isempty(params) || throw(ArgumentError("$name does not take parameters"))
        (:exponential, nothing)
    elseif name in ("gamma", "gamma_edges")
        (:gamma, require_parameter(params, ("shape", "k", "parameter"), "shape"))
    elseif name in ("lognormal", "lognormal_edges")
        (:lognormal, require_parameter(params, ("sigma", "parameter"), "sigma"))
    elseif name in ("pareto", "pareto_edges")
        (:pareto, require_parameter(params, ("alpha", "parameter"), "alpha"))
    elseif name in ("uniform", "uniform_edges")
        (:uniform, require_parameter(params, ("a", "parameter"), "a"))
    elseif name in ("beta", "beta_edges")
        (:beta, require_parameter(params, ("a", "parameter"), "a"))
    elseif name in ("weibull", "weibull_edges")
        (:weibull, require_parameter(params, ("shape", "k", "parameter"), "shape"))
    elseif name in ("inverse_gamma", "inverse_gamma_edges")
        (:inverse_gamma, require_parameter(params, ("alpha", "parameter"), "alpha"))
    elseif name in ("bernoulli", "bernoulli_edges")
        (:bernoulli, require_parameter(params, ("a", "parameter"), "a"))
    else
        (:triangular, require_parameter(params, ("a", "parameter"), "a"))
    end

    if parameter !== nothing
        parameter > 0 || throw(ArgumentError("$name needs a positive parameter"))
        model in (:pareto, :inverse_gamma) && parameter <= 1 &&
            throw(ArgumentError("$name needs a parameter greater than 1"))
        model in (:uniform, :bernoulli, :triangular) && parameter >= 1 &&
            throw(ArgumentError("$name needs 0 < a < 1"))
    end
    return model, parameter
end

function sample_weight(rng::AbstractRNG, model::Symbol, parameter)::Float64
    model === :baseline && return 1.0
    model === :exponential && return rand(rng, Exponential(1.0))
    model === :gamma && return rand(rng, Gamma(parameter, inv(parameter)))
    model === :lognormal && return rand(rng, LogNormal(-0.5 * parameter^2, parameter))
    model === :pareto && return rand(rng, Pareto(parameter, (parameter - 1.0) / parameter))
    model === :uniform && return rand(rng, Uniform(1.0 - parameter, 1.0 + parameter))
    model === :beta && return 2.0 * rand(rng, Beta(parameter, parameter))
    if model === :weibull
        unit = Weibull(parameter, 1.0)
        return rand(rng, Weibull(parameter, inv(mean(unit))))
    end
    model === :inverse_gamma && return rand(rng, InverseGamma(parameter, parameter - 1.0))
    model === :bernoulli && return rand(rng, Bool) ? 1.0 + parameter : 1.0 - parameter
    if model === :triangular
        # Inverse CDF for a symmetric triangular distribution on [1-a, 1+a].
        u = rand(rng)
        return u < 0.5 ? 1.0 - parameter + parameter * sqrt(2u) :
                         1.0 + parameter - parameter * sqrt(2(1.0 - u))
    end
    throw(ArgumentError("unsupported model: $model"))
end

const PARAMETER_NAMES = Dict(
    "gamma" => "shape", "gamma_edges" => "shape",
    "lognormal" => "sigma", "lognormal_edges" => "sigma",
    "pareto" => "alpha", "pareto_edges" => "alpha",
    "uniform" => "a", "uniform_edges" => "a",
    "beta" => "a", "beta_edges" => "a",
    "weibull" => "shape", "weibull_edges" => "shape",
    "inverse_gamma" => "alpha", "inverse_gamma_edges" => "alpha",
    "bernoulli" => "a", "bernoulli_edges" => "a",
    "triangular" => "a", "triangular_edges" => "a",
)

function parse_distribution_argument(text::AbstractString)
    pieces = split(String(text), ':'; limit=2)
    length(pieces) == 1 && return pieces[1], Dict{String,Float64}()
    name = pieces[1]
    haskey(PARAMETER_NAMES, name) || throw(ArgumentError("unknown parameterized distribution: $name"))
    return name, Dict(PARAMETER_NAMES[name] => parse(Float64, pieces[2]))
end
