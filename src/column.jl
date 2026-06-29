function Column(; name, fluid, N_stages, add_reservoir=true, use_ΔT=false)
    
    @named reb = Reboiler(; fluid)    
    stages = [Stage(; name=Symbol(:stage_, j), fluid) for j in 2:N_stages]
    @named cond = CondenserBuffer(; fluid)

    pars = Num[]
    vars = Num[]

    eqs = Equation[]
    push!(eqs, connect(reb.V_out, stages[1].V_in))
    push!(eqs, connect(stages[1].L_out, reb.L_in))
    for j in 1:N_stages-2
        push!(eqs, connect(stages[j].V_out, stages[j+1].V_in))
        push!(eqs, connect(stages[j+1].L_out, stages[j].L_in))
    end
    push!(eqs, connect(stages[N_stages-1].V_out, cond.V_in))
    push!(eqs, connect(cond.L_back, stages[N_stages-1].L_in))

    systems=[cond, reb, stages...]

    if add_reservoir
        @named res = Reservoir(; fluid)
        push!(eqs, connect(cond.L_dist, res.L_in))
        push!(systems, res)
    end

    if use_ΔT
        @variables ΔT(t)
        push!(vars, ΔT)
        push!(eqs, ΔT ~ reb.T - cond.T)
    end

    return ODESystem(eqs, t, vars, pars; name, systems)
end