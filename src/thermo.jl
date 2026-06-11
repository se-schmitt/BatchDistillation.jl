# PH flash
function _ph_flash(fluid, p, h, n)
    # h_safe = isfinite(h) ? h : -30.0 * sum(n) * 1000.0  # rough h_L_sat per mol
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
@register_array_symbolic _ph_flash(fluid, p, h, n) begin 
    size=(1+2*length(model),)
    eltype=Base.promote_eltype(p, h, n)
end

function wrapper_ph_flash(fluid, p, h, n)
    N = length(fluid)
    out = _ph_flash(fluid, p, h, n)
    T = first(out)
    nL = out[2:2+N]
    nV = out[3+N:end]
    return T, nL, nV
end

# bubble temperature
function _bubble_temperature(model, p, x)
    T, _, _, y = bubble_temperature(model, p, x)
    return [T, y...] 
end
@register_array_symbolic _bubble_temperature(fluid, p, x) begin 
    size=(1+length(model),)
    eltype=Base.promote_eltype(p, x)
end

function wrapper_bubble_temperature(fluid, p, x)
    out = _bubble_temperature(fluid, p, x)
    T = first(out)
    y = out[2:end]
    return T, y
end