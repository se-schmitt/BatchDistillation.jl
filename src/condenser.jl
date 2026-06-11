function CondenserBuffer(; name, fluid)
    N = length(fluid)

    @named V_in = MaterialStream()
    @named L_back = MaterialStream()
    @named L_dist = MaterialStream()

    sts = @variables begin
        xB(t)[1:N]      # mole fraction buffer
        xC(t)[1:N]      # mole fraction condenser
        T(t)            # temperature
    end

    pars = @parameters begin
        p 
        nB
        ϵ
    end 

    _T, _ = wrapper_bubble_temperature(fluid, p, xC)

    eqs = [
        T ~ _T,
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
            @. D(xB) ~ V_in.n/nB * (xC - xB)
        )...,
    ]

    return compose(ODESystem(eqs, t, sts, pars; name), [V_in, L_back, L_dist])
end