function Reboiler(; name, fluid)
    N = length(fluid)

    @named L_in  = MaterialStream(; N)
    @named V_out = MaterialStream(; N)

    sts = @variables begin
        nL(t)           # liquid holdup
        T(t)            # temperature
        x(t)[1:N]       # liquid mole fraction
        y(t)[1:N]       # vapor mole fraction
        H(t)            # enthalpy
    end

    pars = @parameters begin
        p
        Q
    end

    _T, _y = wrapper_bubble_temperature(fluid, p, x)

    eqs = Equation[
        # eq. 8
        D(nL) ~ L_in.n - V_out.n,
        # eq. 9
        collect(
            @. D(x) ~ (L_in.n*(L_in.z - x) - V_out.n*(V_out.z - x)) / nL
        )...,
        # eq. 10
        D(H) ~  L_in.n*enthalpy(fluid, p, T, L_in.z; phase=:liquid) -
                V_out.n*enthalpy(fluid, p, T, V_out.z; phase=:vapor) +
                Q,
        # eq. 12
        T ~ _T,
        collect(
            y .~ _y
        )...,
        # eq. 13
        H ~ nL * enthalpy(fluid, p, T, x; phase=:liquid),
        collect(
            V_out.z .~ y
        )...,
        V_out.h ~ enthalpy(fluid, p, T, V_out.z; phase=:vapor),
    ]
    
    return compose(ODESystem(eqs, t, sts, pars; name), [L_in, V_out])
end