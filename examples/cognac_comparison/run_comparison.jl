## setup
using Printf, ADTypes, Suppressor
using MLThermoProperties, Clapeyron, BatchDistillation, GCIdentifier, ChemicalIdentifiers
using ModelingToolkit, OrdinaryDiffEq

path_thermo = joinpath(@__DIR__, "thermo_database")

# numerics
abstol = 1e-6
reltol = 1e-4
dtmax = 50.0
maxiters = 1_000_000
autodiff = AutoFiniteDiff()

# thermodynamic models (comment out to skip)
models = [
    "NRTL"                => (c, pures) -> CompositeModel(c; liquid=NRTL, fluid=pures, liquid_userlocations=["$(path_thermo)/NRTL_unlike.csv"]),
    "modUNIFAC"           => (c, pures) -> CompositeModel(c; liquid=UNIFAC(get_groups_from_name.(c, Ref(UNIFACGroups))), fluid=pures),
    "COGNAC explicit"     => (c, pures) -> CompositeModel(c; liquid=COGNACexp, fluid=pures),
    "HANNA (no ensemble)" => (c, pures) -> CompositeModel(c; liquid=HANNA(c; ensemble=1), fluid=pures),
    "HANNA"               => (c, pures) -> CompositeModel(c; liquid=HANNA, fluid=pures),
]

# mixtures and experiment-specific parameters
components_5comp = ["methanol", "ethanol", "propan-2-ol", "butan-1-ol", "acetone"]
x0_5comp = [0.1919, 0.2321, 0.2242, 0.2414, 0.1104]
mixtures = [
    "3comp" => (; components = ["butan-1-ol", "propan-2-ol", "water"], x0 = [0.429, 0.43, 0.141],
                n_charge = 21.1, tmax = 12000., Q = 230.71, PS = 70000.0, ϵ = 0.30),
    "5comp" => (; components = components_5comp, x0 = x0_5comp,
                n_charge = 1485. / (PR(components_5comp).params.Mw.values' * x0_5comp), tmax = 4000., Q = 350., PS = 50000.0, ϵ = 0.50),
]

# plant-specific parameters (Section 3.4.1)
N_stages = 12
H_col    = 1.5                     # m, packed height -> per-stage height h/(S-1)
ΔT       = 20.0                    # K, head-to-condenser offset (eq. 7)
ΔP       = 93.0                    # Pa, total pressure drop
nB       = 3.3                     # mol, buffer holdup n^B
Tref     = 273.15                  # K
cp_steel = 510.79                  # J/(kg K), structured packing 1.4404
cp_glass = 830.0                   # J/(kg K), borosilicate glass
# masses and correction factors c^j (eq. 8): reboiler vessel vs. column stages
m_steel  = j -> j == 1 ? 0.25  : 0.058
m_glass  = j -> j == 1 ? 1.25  : 0.20
c_corr   = j -> j == 1 ? 1.0   : 2.5
U_mat    = j -> (T -> c_corr(j) * (m_steel(j)*cp_steel + m_glass(j)*cp_glass) * (T - Tref))   # eq. 1 / S26
Q_loss   = j -> j == N_stages ? (T -> -2.0) : (T -> zero(T))
UA_mat   = j -> j == 1 ? 150.0 : 30.0

## run
results = []

@suppress_err for (mix, (; components, x0, n_charge, tmax, Q, PS, ϵ)) in mixtures
    pures = CompositeModel(components,
        liquid=DIPPR105Liquid, liquid_userlocations=["$(path_thermo)/DIPPR105Liquid_like.csv"],
        saturation=DIPPR101Sat,  saturation_userlocations=["$(path_thermo)/DIPPR101Sat_like.csv"],
        gas=TsonopoulosVirial, gas_userlocations=["$(path_thermo)/TsonopoulosVirial_like.csv"],
    )
    P_stage = j -> PS + (N_stages - j) * ΔP / (N_stages - 1)   # eq. 3, j = 1 (reboiler) .. S (head)

    for (name, f_model) in models
        model = f_model(components, pures)

        @named col = Column(; fluid=model, N_stages, use_ΔT=true, stage_h = H_col/(N_stages-1),
                            reb_f_U   = U_mat(1),
                            stage_f_U = j -> U_mat(j),
                            stage_F_Q = j -> Q_loss(j),
                            reb_UA    = UA_mat(1),
                            stage_UA  = j -> UA_mat(j))
        sys = mtkcompile(col)

        init = initialize_total_reflux(sys; fluid=model, N_stages, p=P_stage(1), Q, x_charge=x0, n_charge, f_U = (i, T) -> U_mat(i)(T))

        params = Dict(
            sys.reb.p => P_stage(1), sys.reb.Q => Q,
            sys.cond.p => PS, sys.cond.nB => nB, sys.cond.ϵ => ϵ, sys.ΔT => ΔT,
            (getproperty(sys, Symbol(:stage_, j)).p => P_stage(j) for j in 2:N_stages)...,
        )
        prob = ODEProblem(sys, merge(init.u0, params), (0.0, tmax); guesses = init.guesses, warn_initialize_determined=false)

        reb_nL = ModelingToolkit.getu(sys, sys.reb.nL)
        dryout = DiscreteCallback((u, t, integ) -> reb_nL(integ) < 1.0, terminate!)
        run_solve(prob) = solve(prob, FBDF(; autodiff); callback=dryout, abstol, reltol, dtmax, maxiters)

        run_solve(remake(prob; tspan = (0.0, dtmax)))   # precompilation
        time_solve = @elapsed sol = run_solve(prob)
        time_per_step = time_solve / length(sol.t)
        @printf("%s | %-20s %-10s t_end = %6.0f s, solve = %8.3f s, %.2e s/step\n", mix, name, sol.retcode, sol.t[end], time_solve, time_per_step)
        flush(stdout)

        push!(results, (; mix, name, time_per_step, t = sol.t,
            x_reb = [sol[col.reb.x[j]] for j in eachindex(components)],
            x_buf = [sol[col.cond.xB[j]] for j in eachindex(components)]))
    end
end
@printf("---------------------------------------------------------------\nFinished!")