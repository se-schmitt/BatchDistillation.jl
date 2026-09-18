using Test
using BatchDistillation
using Clapeyron, ModelingToolkit, OrdinaryDiffEq, SciMLBase
using ModelingToolkit: getdefault, unwrap

const FLUID    = UNIFAC(["methanol", "ethanol"])
const P_ATM    = 1.01325e5
const N_STAGES = 5
const CHARGE   = (x_charge = [0.5, 0.5], n_charge = 20.0)
const Q_REB    = 400.0

compile(; N_stages = N_STAGES, kwargs...) =
    mtkcompile(Column(; name = :col, fluid = FLUID, N_stages, kwargs...))

units(sys, N_stages = N_STAGES) =
    [sys.reb, sys.cond, (getproperty(sys, Symbol(:stage_, j)) for j in 2:N_stages)...]

initialize(sys; kwargs...) =
    initialize_total_reflux(sys; fluid = FLUID, N_stages = N_STAGES, p = P_ATM,
                            Q = Q_REB, CHARGE..., kwargs...)

"Operating parameters; `p` now rides on the Column, so it is not repeated here."
opargs(sys) = Dict{Any,Float64}(sys.reb.Q => Q_REB, sys.cond.nB => 0.05, sys.cond.ϵ => 0.2)

const U_WALL = T -> 500.0 * (T - 298.15)

@testset "BatchDistillation" begin

    @testset "column pressure reaches every unit" begin
        for p in (P_ATM, 5.0e4)          # 5e4 differs from the default, so this really checks
            sys = compile(; p)
            @test all(getdefault(unwrap(u.p)) == p for u in units(sys))
        end
        # omitting `p` must keep the pre-existing call signature working
        @test getdefault(unwrap(compile().reb.p)) == 101325.0
    end

    # Each variant must split every unknown into u0 (differential) or guesses
    # (algebraic) with nothing left over, so `fully_determined = true` holds.
    variants = [
        ("plain",   (;),                 (;),                 (;)),
        ("use_ΔT",  (; use_ΔT = true),   (;),                 (; ΔT = 2.0)),
        ("wall/UA", (; reb_UA = 5.0, reb_f_U = U_WALL,
                       stage_UA = j -> 2.0, stage_f_U = j -> U_WALL),
                    (; f_U = (i, T) -> U_WALL(T)),            (;)),
    ]

    @testset "initialization is exactly determined: $name" for (name, colkw, initkw, parkw) in variants
        sys  = compile(; p = P_ATM, colkw...)
        init = initialize(sys; initkw...)

        @test init.retcode == ReturnCode.Success
        @test init.resmax < 1e-8
        @test init.V > 0

        # the u0/guesses split is a partition of the unknowns
        @test length(init.u0) + length(init.guesses) == length(unknowns(sys))
        @test isempty(intersect(keys(init.u0), keys(init.guesses)))

        pars = opargs(sys)
        haskey(parkw, :ΔT) && (pars[sys.ΔT] = parkw.ΔT)
        # fully_determined = true throws if the initialization system is not square
        prob = ODEProblem(sys, merge(Dict(init.u0), pars), (0.0, 60.0);
                          guesses = init.guesses, fully_determined = true)
        sol = solve(prob, FBDF(autodiff = AutoFiniteDiff()); reltol = 1e-6, abstol = 1e-8)
        @test SciMLBase.successful_retcode(sol)
    end

    @testset "README example reproduces" begin
        sys  = compile(; p = P_ATM)
        init = initialize(sys)
        prob = ODEProblem(sys, merge(Dict(init.u0), opargs(sys)), (0.0, 3600.0);
                          guesses = init.guesses, fully_determined = true)
        sol  = solve(prob, FBDF(autodiff = AutoFiniteDiff()); reltol = 1e-6, abstol = 1e-8)

        @test SciMLBase.successful_retcode(sol)
        # the three numbers quoted in README.md
        @test sol[sys.res.nL][end]   ≈ 7.569 rtol = 1e-3
        @test sol[sys.res.x[1]][end] ≈ 0.825 rtol = 1e-3
        @test sol[sys.reb.x[1]][end] ≈ 0.333 rtol = 1e-3

        # methanol is the light boiler: product richer in it than the still it came from
        @test sol[sys.res.x[1]][end] > sol[sys.reb.x[1]][end]

        # Overall mole balance. Internal flows cancel pairwise and the condenser drum
        # holds a constant nB, so reboiler + stages + receiver is conserved exactly.
        inventory(i) = sol[sys.reb.nL][i] + sol[sys.res.nL][i] +
                       sum(sol[getproperty(sys, Symbol(:stage_, j)).nL][i] for j in 2:N_STAGES)
        @test inventory(lastindex(sol.t)) ≈ inventory(1) rtol = 1e-4
        # the charge is the reboiler alone; the stages and drum are filled on top of it
        @test inventory(1) > CHARGE.n_charge
    end

    @testset "column without a reservoir compiles" begin
        # NOTE: initialize_total_reflux seeds sys.res unconditionally, so it currently
        # requires add_reservoir = true; only the flowsheet is checked here.
        @test length(unknowns(compile(; p = P_ATM, add_reservoir = false))) > 0
    end

    @testset "stale utils.jl is gone" begin
        @test !isdefined(BatchDistillation, :build_column)
        @test !isdefined(BatchDistillation, :run_batch)
    end
end
