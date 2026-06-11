module BatchDistillation

using ModelingToolkit, Clapeyron, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: @register_symbolic

const _N_FLOOR = 1e-8
_clip(x) = isfinite(x) ? max(x, _N_FLOOR) : _N_FLOOR

include("thermo.jl")
include("stage.jl")
include("reboiler.jl")
include("condenser.jl")
include("column.jl")
include("connector.jl")
include("utils.jl")

export Column

end # module BatchDistillation
