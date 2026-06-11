# Stage - one equilibrium tray. PH-flash drives equilibrium; soft-drained
# liquid holdup approaches HL_target; vapor leaves fast via τ_V ≪ τ_L.
function Stage(; name, fluid, f_holdup)
    N = length(fluid)

    @named L_in  = MaterialStream()
    @named L_out = MaterialStream()
    @named V_in  = MaterialStream()
    @named V_out = MaterialStream()

    sts = @variables begin
        nL(t)           # liquid holdup
        T(t)            # temperature
        x(t)[1:N]       # liquid mole fraction
        y(t)[1:N]       # vapor mole fraction
        H(t)            # enthalpy
    end

    pars = @parameters begin
        p
    end

    _T, _y = wrapper_bubble_temperature(fluid, p, x)

    eqs = Equation[
        # eq. 8
        D(nL) ~ L_in.n + V_in.n - L_out.n - V_out.n,
        # eq. 9
        collect(
            @. D(x) ~ (L_in.n*(L_in.z - x) + V_in.n*(V_in.z - x) - V_out.n*(V_out.z - x)) / nL
        )...,
        # eq. 10
        D(H) ~  L_in.n*enthalpy(fluid, p, T, L_in.z; phase=:liquid) +
                V_in.n*enthalpy(fluid, p, T, V_in.z; phase=:vapor) - 
                L_out.n*enthalpy(fluid, p, T, L_out.z; phase=:liquid) - 
                V_out.n*enthalpy(fluid, p, T, V_out.z; phase=:vapor),
        # eq. 11
        1 ~ sum(x),
        # eq. 12
        T ~ _T,
        collect(
            y .~ _y
        )...,
        # eq. 13
        H ~ nL * enthalpy(fluid, p, T, x; phase=:liquid),
        # eq. 14
        nL ~ f_holdup(L_in.n)
    ]
    
    return compose(ODESystem(eqs, t, sts, pars; name), [L_in, L_out, V_in, V_out])
end