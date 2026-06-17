function Reservoir(; fluid, name)
    N = length(fluid)

    @named L_in = MaterialStream(; N)

    sts = @variables begin 
        nL(t)
        x(t)[1:N]
    end

    eqs = Equation[
        D(nL) ~ L_in.n,
        collect(
            @. D(x) ~ (L_in.n*(L_in.z - x)) / nL
        )...,
    ]

    return compose(ODESystem(eqs, t, sts, []; name), [L_in])
end