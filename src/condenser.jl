function CondenserBuffer(; name, fluid, p_stage=101325.0, use_ΔT=false)
    N = length(fluid)

    @named V_in = MaterialStream(; N)
    @named L_back = MaterialStream(; N)
    @named L_dist = MaterialStream(; N)

    vars = @variables begin
        xB(t)[1:N]      # mole fraction buffer
        xC(t)[1:N]      # mole fraction condenser
        T(t)            # temperature
    end

    pars = @parameters begin
        p = p_stage
        nB
        ϵ
    end

    eqs = [
        collect(
            L_dist.z .~ xB
        )...,
        collect(
            L_back.z .~ xB
        )...,
        collect(
            V_in.z .~ xC
        )...,
        L_dist.n ~ ϵ*V_in.n,
        L_back.n ~ (1-ϵ)*V_in.n,
        collect(
            D(xB) .~ (V_in.n/nB) .* (xC .- xB)
        )...,
        L_back.h ~ h_molar(fluid, p, T, L_back.z; phase=:liquid),
        L_dist.h ~ h_molar(fluid, p, T, L_dist.z; phase=:liquid),
    ]

    if !use_ΔT
        _T, _ = wrapper_bubble_temperature(fluid, p, xC)
        push!(eqs, T ~ _T)
    end

    return compose(ODESystem(eqs, t, vars, pars; name), [V_in, L_back, L_dist])
end