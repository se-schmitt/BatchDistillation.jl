function Column(; name, )
    @named cond = Condenser(; R)
    @named reb  = Reboiler(;  Q_reb, τ_V)
    stages = [Stage(; name=Symbol(:stage_, j), HL_target, τ_L, τ_V) for j in 1:N]

    conns = Equation[]
    push!(conns, connect(cond.reflux_out,    stages[1].liq_in))
    push!(conns, connect(stages[1].vap_out,  cond.vap_in))
    for j in 1:N-1
        push!(conns, connect(stages[j].liq_out,   stages[j+1].liq_in))
        push!(conns, connect(stages[j+1].vap_out, stages[j].vap_in))
    end
    push!(conns, connect(stages[N].liq_out, reb.liq_in))
    push!(conns, connect(reb.vap_out,       stages[N].vap_in))

    @named col = ODESystem(conns, t; systems=[cond, reb, stages...])
end