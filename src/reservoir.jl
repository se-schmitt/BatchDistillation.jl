function Reservoir(; fluid, name)
    N = length(fluid)

    @named L_in = MaterialStream(; N)

    sts = @variables begin 
        nL(t)
        x(t)[1:N]
    end

    eqs = Equation[
        D(nL) ~ L_in.n,
        [D(x[i]) ~ (L_in.n*(L_in.z[i] - x[i])) / nL for i in 1:N-1]...,
        1 ~ sum(collect(x)),
    ]

    return compose(ODESystem(eqs, t, sts, []; name), [L_in])
end