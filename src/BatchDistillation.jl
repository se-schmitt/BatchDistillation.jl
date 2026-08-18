module BatchDistillation

using ModelingToolkit, Clapeyron, OrdinaryDiffEq, NonlinearSolve
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: @register_symbolic

const _N_FLOOR = 1e-8
_clip(x) = isfinite(x) ? (_N_FLOOR + 0.5*(x + sqrt(x*x + _N_FLOOR*_N_FLOOR))) : _N_FLOOR

# Thermo
include("thermo.jl")

# Components
include("stage.jl")
include("reboiler.jl")
include("condenser.jl")
include("column.jl")
include("connector.jl")
include("reservoir.jl")

# Utils
include("utils.jl")
include("initialization.jl")

export Column, initialize_total_reflux

end # module BatchDistillation
