function Column(; name, fluid, N_stages, add_reservoir=true, use_ΔT=false,
                stage_h = 1.0,
                reb_f_U = xT->zero(xT),
                stage_f_U = j -> (xT->zero(xT)),
                stage_F_Q = j -> (xT->zero(xT)))

    @named reb = Reboiler(; fluid, f_U_stage=reb_f_U)
    stages = [Stage(; name=Symbol(:stage_, j), fluid, h_stage=stage_h,
                    f_U_stage=stage_f_U(j), F_Q_stage=stage_F_Q(j)) for j in 2:N_stages]
    @named cond = CondenserBuffer(; fluid, use_ΔT)

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
        @parameters ΔT
        push!(pars, ΔT)
        push!(eqs, cond.T ~ stages[N_stages-1].T - ΔT)
    end

    return ODESystem(eqs, t, vars, pars; name, systems)
end