# Assemble the full column: condenser on top, N stages, reboiler at bottom.
function build_column(; N=3, R=2.0, Q_reb=500.0,
                       HL_target=1.0, τ_L=1.0, τ_V=1e-3)
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
    return col, cond, reb, stages
end

# Define initial conditions
function initial_conditions(cond, reb, stages;
        n_init=0.1, n0_reservoir=10.0, h_init_offset_J_per_mol=2000.0)
    z_eq = [1/3, 1/3, 1/3]
    T_bub_init = bubble_temperature(model, p_set, z_eq)[1]
    h_L_init = enthalpy(model, p_set, T_bub_init, z_eq; phase=:liquid)
    h_init = h_L_init + h_init_offset_J_per_mol
    u0 = Pair[]
    push!(u0, reb.n1 => n0_reservoir/3)
    push!(u0, reb.n2 => n0_reservoir/3)
    push!(u0, reb.n3 => n0_reservoir/3)
    push!(u0, reb.H_stage => n0_reservoir * h_init)
    for s in stages
        push!(u0, s.n1 => n_init/3)
        push!(u0, s.n2 => n_init/3)
        push!(u0, s.n3 => n_init/3)
        push!(u0, s.H_stage => n_init * h_init)
    end
    push!(u0, cond.n_dist1 => 0.0)
    push!(u0, cond.n_dist2 => 0.0)
    push!(u0, cond.n_dist3 => 0.0)
    return u0
end

# End-to-end run: assemble, simplify, solve. Returns (sol, cond, reb, stages).
function run_batch(; N=3, R=2.0, Q_reb=500.0, tend=600.0,
                    HL_target=1.0, τ_L=1.0, τ_V=0.1,
                    n_init=0.1, n0_reservoir=10.0, h_init_offset_J_per_mol=2000.0,
                    reltol=1e-4, abstol=1e-6, dtmax=5.0)
    col, cond, reb, stages = build_column(; N, R, Q_reb, HL_target, τ_L, τ_V)
    sys = structural_simplify(col)
    u0 = initial_conditions(cond, reb, stages; n_init, n0_reservoir, h_init_offset_J_per_mol)
    prob = ODEProblem(sys, u0, (0.0, tend))
    sol = solve(prob, FBDF(autodiff=false); reltol, abstol, dtmax)
    return sol, cond, reb, stages
end