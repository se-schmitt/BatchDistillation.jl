# Condenser - total condenser at the bubble point of the condensate. Splits
# condensate into reflux (L = R·D) and distillate (D). Accumulates n_dist_i.
function Condenser(; name)
    N = length(fluid)

    @named V_in = MaterialStream()
    @named L_back = MaterialStream()
    @named L_dist = MaterialStream()

    sts = @variables begin
        xB[1:N](t)      # mole fraction buffer
        xC[1:N](t)      # mole fraction condenser
        T(t)            # temperature
    end
    pars = @parameters begin
        p 
        nB
        ϵ
    end 

    _T, _ = wrapper_bubble_temperature(fluid, p, V_in.z)

    eqs = [
        T ~ _T,
        collect(
            L_dist.z .~ xB
        )...,
        collect(
            L_back.z .~ xB
        )...,
        L_dist.n ~ ϵ*V_in.n,
        L_back.n ~ (1-ϵ)*V_in.n,
        collect(
            @. D(xB) ~ V_in.n/nB * (xC - xB)
        )...,
    ]

    return compose(ODESystem(eqs, t, sts, pars; name), [vap_in, reflux_out])
end