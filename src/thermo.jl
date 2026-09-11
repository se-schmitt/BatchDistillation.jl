using ModelingToolkit.Symbolics: @register_derivative
const ForwardDiff = Clapeyron.ForwardDiff

# clamp to a positive, normalized composition so Clapeyron never sees infeasible iterates
_safez(xs::Vararg) = (v = map(_clip, collect(xs)); v ./ sum(v))

_pure_models(m) = hasproperty(m, :puremodel) ? m.puremodel.pure :
                  hasproperty(m, :vle)       ? _pure_models(m.vle) :
                  hasproperty(m, :fluid)     ? _pure_models(m.fluid) : nothing

const _TSAT_CACHE = Dict{Tuple{UInt,Int,Float64},Float64}()

function _T0_guess(fluid, p, z)
    pure = _pure_models(fluid)
    isnothing(pure) && return nothing
    pv = Float64(ForwardDiff.value(p))
    id = objectid(pure)
    T0 = 0.0
    for i in eachindex(z)
        Ts = get!(_TSAT_CACHE, (id, i, pv)) do
            Float64(Clapeyron.saturation_temperature(pure[i], pv)[1])
        end
        isfinite(Ts) || return nothing
        T0 += Float64(ForwardDiff.value(z[i])) * Ts
    end
    return T0
end

# bubble_temperature with initial guess
function _bubbleT(fluid, p, z)
    T0 = _T0_guess(fluid, p, z)
    return isnothing(T0) ? bubble_temperature(fluid, p, z) :
                           bubble_temperature(fluid, p, z; T0 = T0)
end

function _ph_flash(fluid, p, h, n)
    z = n ./ sum(n)
    res = ph_flash(fluid, p, h, n)
    T = Clapeyron.temperature(res)
    if length(res.compositions) == 1
        nL = res.compositions[1] .* res.fractions[1]
        return [T, nL, zeros(length(fluid))...]
    else
        nL = res.compositions[1] .* res.fractions[1]
        nV = res.compositions[2] .* res.fractions[2]
        return [T, nL..., nV...]
    end
end
@register_array_symbolic _ph_flash(fluid, p, h, n::AbstractVector) begin
    size=(1+2*length(model),)
    eltype=eltype(n)
    ndims=1
end

function wrapper_ph_flash(fluid, p, h, n)
    N = length(fluid)
    out = _ph_flash(fluid, p, h, n)
    T = first(out)
    nL = out[2:2+N]
    nV = out[3+N:end]
    return T, nL, nV
end

# Scalarized, ForwardDiff-differentiable enthalpy and bubble-temperature wrappers
# (composition as scalar components so MTK index reduction can differentiate them).
for N in 2:5
    xs = [Symbol(:x, i) for i in 1:N]
    bT = Symbol(:_bT_, N)

    for (suf, ph) in ((:L, QuoteNode(:liquid)), (:V, QuoteNode(:vapor)))
        enth = Symbol(:_enth, suf, :_, N)
        @eval begin
            $enth(fluid, p, T, $(xs...)) = enthalpy(fluid, p, T, _safez($(xs...)); phase=$ph)
            @register_symbolic $enth(fluid::Clapeyron.EoSModel, p, T, $(xs...))
            $(Symbol(enth, :_dT))(fluid, p, T, $(xs...)) =
                ForwardDiff.derivative(θ -> $enth(fluid, p, θ, $(xs...)), T)
            @register_symbolic $(Symbol(enth, :_dT))(fluid::Clapeyron.EoSModel, p, T, $(xs...))
            @register_derivative $enth(fluid, p, T, $(xs...)) 3 $(Symbol(enth, :_dT))(fluid, p, T, $(xs...))
        end
        for i in 1:N
            d = Symbol(enth, :_dx, i)
            @eval begin
                $d(fluid, p, T, $(xs...)) =
                    ForwardDiff.derivative(ξ -> $enth(fluid, p, T, $(xs[1:i-1]...), ξ, $(xs[i+1:N]...)), $(xs[i]))
                @register_symbolic $d(fluid::Clapeyron.EoSModel, p, T, $(xs...))
                @register_derivative $enth(fluid, p, T, $(xs...)) $(3 + i) $d(fluid, p, T, $(xs...))
            end
        end
    end

    @eval begin
        $bT(fluid, p, $(xs...)) = _bubbleT(fluid, p, _safez($(xs...)))[1]
        @register_symbolic $bT(fluid::Clapeyron.EoSModel, p, $(xs...))
    end
    for i in 1:N
        d = Symbol(bT, :_dx, i)
        @eval begin
            $d(fluid, p, $(xs...)) =
                ForwardDiff.derivative(ξ -> $bT(fluid, p, $(xs[1:i-1]...), ξ, $(xs[i+1:N]...)), $(xs[i]))
            @register_symbolic $d(fluid::Clapeyron.EoSModel, p, $(xs...))
            @register_derivative $bT(fluid, p, $(xs...)) $(2 + i) $d(fluid, p, $(xs...))
        end
    end

    for k in 1:N
        bY = Symbol(:_bY_, N, :_, k)
        @eval begin
            $bY(fluid, p, $(xs...)) = _bubbleT(fluid, p, _safez($(xs...)))[4][$k]
            @register_symbolic $bY(fluid::Clapeyron.EoSModel, p, $(xs...))
        end
        for i in 1:N
            d = Symbol(bY, :_dx, i)
            @eval begin
                $d(fluid, p, $(xs...)) =
                    ForwardDiff.derivative(ξ -> $bY(fluid, p, $(xs[1:i-1]...), ξ, $(xs[i+1:N]...)), $(xs[i]))
                @register_symbolic $d(fluid::Clapeyron.EoSModel, p, $(xs...))
                @register_derivative $bY(fluid, p, $(xs...)) $(2 + i) $d(fluid, p, $(xs...))
            end
        end
    end
end

function h_molar(fluid, p, T, z; phase)
    suf = phase === :liquid ? :L : :V
    f = getfield(@__MODULE__, Symbol(:_enth, suf, :_, length(z)))
    f(fluid, p, T, collect(z)...)
end

function wrapper_bubble_temperature(fluid, p, x)
    N = length(x)
    xs = collect(x)
    T = getfield(@__MODULE__, Symbol(:_bT_, N))(fluid, p, xs...)
    y = [getfield(@__MODULE__, Symbol(:_bY_, N, :_, k))(fluid, p, xs...) for k in 1:N]
    return T, y
end
