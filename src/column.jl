function Column(; name, fluid, N_stages)
    
    @named reb = Reboiler(; fluid)    
    stages = [Stage(; name=Symbol(:stage_, j), fluid) for j in 2:N_stages]
    @named cond = CondenserBuffer(; fluid)

    conns = Equation[]
    push!(conns, connect(reb.V_out, stages[1].V_in))
    push!(conns, connect(stages[1].L_out, reb.L_in))
    for j in 1:N_stages-2
        push!(conns, connect(stages[j].V_out, stages[j+1].V_in))
        push!(conns, connect(stages[j+1].L_out, stages[j].L_in))
    end
    push!(conns, connect(stages[N_stages-1].V_out, cond.V_in))
    push!(conns, connect(cond.L_back, stages[N_stages-1].L_in))

    return ODESystem(conns, t; name, systems=[cond, reb, stages...])
end