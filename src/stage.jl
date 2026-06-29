function Stage(; name, fluid, f_U_stage=xT->zero(xT))
    N = length(fluid)

    @named L_in  = MaterialStream(; N)
    @named L_out = MaterialStream(; N)
    @named V_in  = MaterialStream(; N)
    @named V_out = MaterialStream(; N)

    vars = @variables begin
        nL(t)           # liquid holdup
        T(t)            # temperature
        x(t)[1:N]       # liquid mole fraction
        y(t)[1:N]       # vapor mole fraction
        U(t)            # enthalpy
    end

    pars = @parameters begin
        p
        n_ref = 0.042
        B_ref = 5/3600
        h = 1.0
    end

    _T, _y = wrapper_bubble_temperature(fluid, p, x)

    eqs = Equation[
        # eq. 8
        D(nL) ~ L_in.n + V_in.n - L_out.n - V_out.n,
        # eq. 9 (S24)
        [D(x[i]) ~ (L_in.n*(L_in.z[i] - x[i]) + V_in.n*(V_in.z[i] - x[i]) - V_out.n*(V_out.z[i] - x[i])) / nL for i in 1:N-1]...,
        1 ~ sum(collect(x)),
        # eq. 10
        D(U) ~  L_in.n*L_in.h + V_in.n*V_in.h - L_out.n*L_out.h - V_out.n*V_out.h,
        # eq. 12
        T ~ _T,
        collect(
            y .~ _y
        )...,
        # eq. 13
        U ~ nL * h_molar(fluid, p, T, x; phase=:liquid) + f_U_stage(T),
        collect(
            L_out.z .~ x
        )...,
        collect(
            V_out.z .~ y
        )...,
        L_out.h ~ h_molar(fluid, p, T, L_out.z; phase=:liquid),
        V_out.h ~ h_molar(fluid, p, T, V_out.z; phase=:vapor),
        # eq. 14 (S27)
        nL ~ n_ref / B_ref * h * L_out.n
    ]
    
    return compose(ODESystem(eqs, t, vars, pars; name), [L_in, L_out, V_in, V_out])
end