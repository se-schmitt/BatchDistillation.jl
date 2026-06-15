# Total-reflux (ε=0) steady-state init -> (; sys, u0): consistent start (free states +
# flow guesses, no determined T/H) for a transient. Σx=1 is re-imposed as extra
# least-squares residuals to pin the per-stage composition scale dropped for DOF balance.
function initialize_total_reflux(; fluid, N_stages, p, Q, nB, x_charge,
                                 n_charge = 20.0, nL_stage = 0.04, n_res = 1e-3,
                                 V_guess = Q / 4e4, abstol = 1e-10, maxiters = 1000)
    @named col = Column(; fluid, N_stages)
    sys = mtkcompile(col)
    N = length(fluid)
    setc!(d, v, val) = for i in 1:N; d[v[i]] = val[i]; end

    gp = Dict{Any,Any}(sys.reb.p => p, sys.reb.Q => Q,
                       sys.cond.p => p, sys.cond.nB => nB, sys.cond.ϵ => 0.0)
    gp[sys.reb.nL] = n_charge; setc!(gp, sys.reb.x, x_charge)
    setc!(gp, sys.cond.xB, x_charge)
    gp[sys.res.nL] = n_res; setc!(gp, sys.res.x, x_charge)
    stages = [getproperty(sys, Symbol(:stage_, j)) for j in 2:N_stages]
    for st in stages
        gp[st.p] = p; gp[st.nL] = nL_stage; setc!(gp, st.x, x_charge)
    end

    odeprob = ODEProblem(sys, gp, (0.0, 1.0))
    U = unknowns(sys); uw = ModelingToolkit.unwrap.(U)
    idx(v) = findfirst(isequal(ModelingToolkit.unwrap(v)), uw)
    isflow(v) = occursin(r"₊(V_out|V_in|L_out|L_in|L_back|L_dist)₊n", string(v))
    u_start = copy(odeprob.u0)
    for i in eachindex(U); isflow(U[i]) && (u_start[i] = V_guess); end

    compvecs = [sys.reb.x, sys.cond.xB, sys.res.x, (st.x for st in stages)...]
    groups = [[idx(v[i]) for i in 1:N] for v in compvecs]
    m = length(U); nres = m + length(groups)
    function resid!(out, u, par)
        odeprob.f(view(out, 1:m), u, par, 0.0)
        for (k, g) in enumerate(groups); out[m+k] = sum(u[i] for i in g) - 1; end
        out
    end
    nlf = NonlinearFunction(resid!; resid_prototype = zeros(nres))
    nlsol = solve(NonlinearLeastSquaresProblem(nlf, u_start, odeprob.p),
                  LevenbergMarquardt(); abstol, maxiters)

    isbound(v) = occursin("₊T(", string(v)) || occursin("₊H(", string(v))
    u0 = Dict(U[i] => nlsol.u[i] for i in eachindex(U) if !isbound(U[i]))
    return (; sys, u0, retcode = nlsol.retcode, resmax = maximum(abs, nlsol.resid))
end
