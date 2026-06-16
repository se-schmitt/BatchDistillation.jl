function Reservoir(; fluid, name)
    N = length(fluid)

    @named L_in = MaterialStream(; N)

    sts = @variables begin 
        nL(t)
        x(t)[1:N]
    end

    dxdt = collect(@. (L_in.n*(L_in.z - x)) / nL)

    eqs = Equation[
        D(nL) ~ L_in.n,
        [D(x[i]) ~ dxdt[i] for i in 1:N-1]...,
        1 ~ sum(x),
    ]

    return compose(ODESystem(eqs, t, sts, []; name), [L_in])
end