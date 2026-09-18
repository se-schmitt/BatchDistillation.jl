## setup
using Printf, Statistics, ADTypes, DataFrames, CSV
using Clapeyron, BatchDistillation, GCIdentifier, ChemicalIdentifiers
using ModelingToolkit, SciMLBase, OrdinaryDiffEq
using CairoMakie

# With no argument the outputs land next to this script (reproducing the committed
# figure); pass a name to write into a subfolder instead.
run_name     = isempty(ARGS) ? "" : ARGS[1]
path_results = isempty(run_name) ? (@__DIR__) : joinpath(@__DIR__, run_name)
mkpath(path_results)

# numerics
abstol = 1e-6
reltol = 1e-4
dtmax = 50.0
maxiters = 1_000_000
autodiff = AutoFiniteDiff()

# input
tmax = 13000.
n_dry  = 0.5                       # mol, reboiler holdup at which the batch is stopped
x_pure = 1e-4                      # stop once the reboiler reaches 1 - x_pure in any component

# system
components = ["butan-1-ol", "propan-2-ol", "water"]

# thermodynamic models
groups_modunifac = get_groups_from_name.(components, Ref(UNIFACGroups))
model = UNIFAC(groups_modunifac; puremodel = SRK)

# g/mol; the activity/ML models carry no Mw
Mw = PR(components).params.Mw.values

# plant-specific parameters (Section 3.4.1)
N_stages = 12
H_col    = 1.5                     # m, packed height -> per-stage height h/(S-1)
ΔT       = 20.0                    # K, head-to-condenser offset (eq. 7)
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

# experiment-specific parameters (calibration experiment, Table 1, ID 1)
Q    = 230.71                      # W, reboiler heat duty (already net of the reboiler heat loss, eq. 9)
PS   = 70000.0                     # Pa, head pressure P^S
ΔP   = 93.0                        # Pa, total pressure drop
ϵ        = 0.30                        # efflux ratio
x0   = [0.429, 0.43, 0.141]        # reboiler charge x^1(0)  (butanol, 2-propanol, water)
nB   = 3.3                         # mol, buffer holdup n^B (SI S1.2: 0.075 kg vessel + 95 mL tubing)
n_charge = 21.1                    # mol, reboiler charge so that n^app(0) ≈ 24.62 (Table 1)

P_stage(j) = PS + (N_stages - j) * ΔP / (N_stages - 1)   # eq. 3, j = 1 (reboiler) .. S (head)

## run
# create column
@named col = Column(; fluid=model, N_stages, use_ΔT=true, stage_h = H_col/(N_stages-1),
                reb_f_U   = U_mat(1),
                stage_f_U = j -> U_mat(j),
                stage_F_Q = j -> Q_loss(j),
                reb_UA    = UA_mat(1),
                stage_UA  = j -> UA_mat(j))
sys = mtkcompile(col)

# parameter map: duty, buffer, efflux ratio, head offset and the pressure profile
params = Dict(
	sys.reb.p => P_stage(1), sys.reb.Q => Q,
	sys.cond.p => PS, sys.cond.nB => nB, sys.cond.ϵ => ϵ, sys.ΔT => ΔT,
	(getproperty(sys, Symbol(:stage_, j)).p => P_stage(j) for j in 2:N_stages)...,
)

# initialization
init = initialize_total_reflux(sys; fluid=model, N_stages, p=P_stage(1), Q, x_charge=x0, n_charge, f_U = (i, T) -> U_mat(i)(T))

prob = ODEProblem(sys, merge(init.u0, params), (0.0, tmax); guesses = init.guesses, warn_initialize_determined=false)

# callbacks/events
reb_nL = ModelingToolkit.getu(sys, sys.reb.nL)
reb_x  = ModelingToolkit.getu(sys, collect(sys.reb.x))

# stop the batch once the reboiler runs dry, or once it is essentially a pure
# component — past that the composition leaves the simplex and the solve fails
# save_positions = (true, false): without it the terminating save appends a
# duplicate point carrying the initial state.
_dryout = DiscreteCallback(
    (u, t, integ) -> reb_nL(integ) < n_dry || maximum(reb_x(integ)) > 1 - x_pure,
    terminate!; save_positions = (true, false))

# solve
sol = solve(prob, FBDF(; autodiff); callback=_dryout, abstol, reltol, dtmax, maxiters)

@info "solve finished" sol.retcode t_end_s=sol.t[end] n_reb_end=sol[sys.reb.nL][end]

## Visualization of results

# --- pull the requested series -------------------------------------------------
t_min  = sol.t ./ 60
T_reb  = sol[sys.reb.T]
T_cond = sol[sys.cond.T]
x_reb  = [sol[sys.reb.x[i]]   for i in eachindex(components)]   # reboiler liquid
x_cond = [sol[sys.cond.xB[i]] for i in eachindex(components)]   # condenser/buffer liquid

stopped_early = sol.t[end] < tmax - 1.0
dried_out     = sol.retcode == SciMLBase.ReturnCode.Terminated   # the _dryout callback fired

# table view (also the relief for the low-contrast series colour)
df = DataFrame(t_s = sol.t, t_min = t_min, T_reb_K = T_reb, T_cond_K = T_cond)
for (i, c) in pairs(components)
    df[!, Symbol("x_reb_",  replace(c, "-" => "_"))] = x_reb[i]
    df[!, Symbol("x_cond_", replace(c, "-" => "_"))] = x_cond[i]
end
CSV.write(joinpath(path_results, "profiles.csv"), df)

# --- figure --------------------------------------------------------------------
# categorical slots 1-3, validated for all pairs in light mode
series_colors = ["#2a78d6", "#eb6834", "#1baf7a"]
ink_primary, ink_muted = "#0b0b0b", "#52514e"

set_theme!(Theme(
    fontsize = 14,
    Axis = (; backgroundcolor = :transparent,
              xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
              xgridwidth = 1, ygridwidth = 1,
              leftspinevisible = true, bottomspinevisible = true,
              topspinevisible = false, rightspinevisible = false,
              spinewidth = 1, xtickcolor = ink_muted, ytickcolor = ink_muted,
              xlabelcolor = ink_primary, ylabelcolor = ink_primary,
              titlecolor = ink_primary, titlealign = :left, titlesize = 15),
))

fig = Figure(size = (980, 640), backgroundcolor = :white)

ax_Treb  = Axis(fig[1, 1]; title = "Reboiler",  ylabel = "T / K")
ax_Tcond = Axis(fig[1, 2]; title = "Condenser", ylabel = "T / K")
ax_xreb  = Axis(fig[2, 1]; ylabel = "x / mol mol⁻¹", xlabel = "t / min")
ax_xcond = Axis(fig[2, 2]; ylabel = "x / mol mol⁻¹", xlabel = "t / min")

lines!(ax_Treb,  t_min, T_reb;  color = ink_primary, linewidth = 2)
lines!(ax_Tcond, t_min, T_cond; color = ink_primary, linewidth = 2)

for (i, c) in pairs(components)
    lines!(ax_xreb,  t_min, x_reb[i];  color = series_colors[i], linewidth = 2, label = c)
    lines!(ax_xcond, t_min, x_cond[i]; color = series_colors[i], linewidth = 2, label = c)
end

# direct labels at the line ends (relief for the low-contrast slot).
# nudge apart any labels whose series end too close together.
function spread(vals; mingap)
    p, out = sortperm(vals), copy(vals)
    for k in 2:length(p)
        out[p[k]] = max(out[p[k]], out[p[k-1]] + mingap)
    end
    out
end
for (ax, xs) in ((ax_xreb, x_reb), (ax_xcond, x_cond))
    ylab = spread([x[end] for x in xs]; mingap = 0.085)
    for (i, c) in pairs(components)
        text!(ax, t_min[end], ylab[i]; text = " " * c, color = ink_muted,
              align = (:left, :center), fontsize = 11)
    end
end

# mark where the run ended early, and say why (label sits in the right-hand padding)
if stopped_early
    why = if !dried_out
        "stopped: $(sol.retcode)"
    elseif sol[sys.reb.nL][end] < n_dry
        @sprintf("dry-out, n¹ = %.1f mol", n_dry)
    else
        @sprintf("reboiler at %g %% purity", 100 * (1 - x_pure))
    end
    for ax in (ax_Treb, ax_Tcond, ax_xreb, ax_xcond)
        vlines!(ax, t_min[end]; color = (ink_muted, 0.45), linewidth = 1, linestyle = :dash)
    end
    text!(ax_Treb, t_min[end], mean(extrema(T_reb)); text = why, color = ink_muted,
          align = (:center, :bottom), fontsize = 11, rotation = π/2, offset = (4, 0))
end

for ax in (ax_Treb, ax_Tcond, ax_xreb, ax_xcond)
    xlims!(ax, 0, t_min[end] * 1.22)          # room for the direct labels
end
for ax in (ax_Treb, ax_Tcond)
    hidexdecorations!(ax; grid = false)
end
linkxaxes!(ax_Treb, ax_Tcond, ax_xreb, ax_xcond)
linkyaxes!(ax_xreb, ax_xcond)

Legend(fig[3, 1:2], ax_xreb; orientation = :horizontal, framevisible = false,
       labelcolor = ink_primary, padding = (0, 0, 0, 0))

Label(fig[0, 1:2], "Batch distillation: butan-1-ol / propan-2-ol / water";
      fontsize = 16, color = ink_primary, halign = :left, padding = (4, 0, 0, 0))

rowgap!(fig.layout, 8)
colgap!(fig.layout, 26)

save(joinpath(path_results, "T_x_reboiler_condenser.png"), fig; px_per_unit = 2)
@info "results written" path_results
