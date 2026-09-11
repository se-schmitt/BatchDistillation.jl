#TODO replace core by adiabatic column component 
function initialize_total_reflux(sys; fluid, N_stages, p, Q, x_charge,
                                 n_charge = 20.0, n_res = 1e-3,
                                 f_U = (i, T) -> zero(float(T)),
                                 abstol = 1e-10, maxiters = 1000)
    N = length(fluid)
    M = N_stages                      # liquid units: 1 = reboiler, 2..M = stages

    _vapor(x) = bubble_temperature(fluid, p, _safez(x...))[4]

    # reduced replacement system: M·N composition unknowns (flows cancel at total reflux)
    function resid!(F, X, _)
        xs = [view(X, (i-1)*N+1 : i*N) for i in 1:M]
        F[1:N] .= xs[1] .- x_charge                  # pin reboiler to charge
        for i in 2:M
            y = _vapor(xs[i-1])                       # passing vapor from the unit below
            b = (i-1)*N
            for j in 1:N-1
                F[b+j] = xs[i][j] - y[j]              # x_i = y_{i-1} (drop one component)
            end
            F[b+N] = sum(xs[i]) - 1                   # Σx = 1
        end
        F
    end
    X0 = Vector{Float64}(undef, M*N)
    X0[1:N] .= float.(x_charge)
    for i in 2:M
        X0[(i-1)*N+1 : i*N] .= _vapor(view(X0, (i-2)*N+1 : (i-1)*N))
    end
    nlsol = solve(NonlinearProblem(NonlinearFunction(resid!), X0), NewtonRaphson(); abstol, maxiters)

    xcol = [nlsol.u[(i-1)*N+1 : i*N] for i in 1:M]    # liquid comp per unit
    Tcol = [bubble_temperature(fluid, p, _safez(xc...))[1] for xc in xcol]
    ycol = [_vapor(xc) for xc in xcol]                # vapor comp per unit
    xC   = ycol[M]                                    # vapor to condenser = top vapor
    
    hVout = h_molar(fluid, p, Tcol[1], ycol[1]; phase = :vapor)
    hLin  = M ≥ 2 ? h_molar(fluid, p, Tcol[2], xcol[2]; phase = :liquid) :
                    h_molar(fluid, p, Tcol[1], xcol[1]; phase = :liquid)
    V = Q / (hVout - hLin)

    setv!(d, v, val) = (d[ModelingToolkit.unwrap(v)] = val)
    getp(v) = ModelingToolkit.getdefault(ModelingToolkit.unwrap(v))
    vals = Dict{Any,Float64}()
    units = [sys.reb, (getproperty(sys, Symbol(:stage_, j)) for j in 2:N_stages)...]
    for (i, u) in enumerate(units)
        # reboiler holdup is the pinned charge; stage holdups follow nL = k·L_out.n (=k·V at total reflux)
        nL = i == 1 ? n_charge : getp(u.n_ref) / getp(u.B_ref) * getp(u.h) * V
        for k in 1:N
            setv!(vals, u.x[k], xcol[i][k])
            hasproperty(u, :L_out) && setv!(vals, u.L_out.z[k], xcol[i][k])
            hasproperty(u, :V_out) && setv!(vals, u.V_out.z[k], ycol[i][k])
        end
        setv!(vals, u.nL, nL)
        setv!(vals, u.T, Tcol[i])
        wall = hasproperty(u, :T_wall)
        setv!(vals, u.U, nL * h_molar(fluid, p, Tcol[i], xcol[i]; phase = :liquid) +
                         (wall ? zero(Tcol[i]) : f_U(i, Tcol[i])))
        if wall
            setv!(vals, u.T_wall, Tcol[i])
            setv!(vals, u.U_wall, f_U(i, Tcol[i]))
        end
        if hasproperty(u, :L_in)
            zin = i < M ? xcol[i+1] : xC
            for k in 1:N; setv!(vals, u.L_in.z[k], zin[k]); end
        end
        if hasproperty(u, :V_in)
            for k in 1:N; setv!(vals, u.V_in.z[k], ycol[i-1][k]); end
        end
    end
    for k in 1:N
        setv!(vals, sys.cond.xB[k], xC[k])
        setv!(vals, sys.res.x[k], x_charge[k])
    end
    setv!(vals, sys.res.nL, n_res)

    isflow(s)  = occursin(r"₊(V_out|V_in|L_out|L_in|L_back|L_dist)₊n", s)
    isderiv(s) = occursin("ˍt", s)                  # dummy-derivative state (index reduction)
    isalg(s)   = isflow(s) || isderiv(s) || occursin(r"₊[TU]\(", s) || occursin(r"stage_\d+₊nL", s)
    u0      = Dict{Any,Float64}()
    guesses = Dict{Any,Float64}()
    for u in unknowns(sys)
        uw = ModelingToolkit.unwrap(u); s = string(u)
        val = haskey(vals, uw) ? vals[uw] : isflow(s) ? V : 0.0
        (isalg(s) ? guesses : u0)[u] = val
    end

    return (; sys, u0, guesses, V, retcode = nlsol.retcode, resmax = maximum(abs, nlsol.resid))
end
