# COGNAC comparison

This script accompanies the paper *Learning atomic interactions for the thermodynamics of liquid mixtures* by Hoffmann et al., [arXiv:2610.36440](https://arxiv.org/abs/2610.36440) (2026). The script simulates one batch distillation model with different activity coefficient models and measures the solve time of each simulation.

The simulation closely follows the batch distillation model described in detail by Werner et al., *Automated batch distillation process simulation for a large hybrid dataset for deep anomaly detection*, Comput. Chem. Eng. 215 (2026) 109827 ([doi:10.1016/j.compchemeng.2026.109827](https://doi.org/10.1016/j.compchemeng.2026.109827)).

[`run_comparison.jl`](run_comparison.jl) simulates a 12-unit packed batch column (reboiler + 11 stages, wall heat capacities, heat loss at the top, linear pressure drop) for two mixtures:

- three component mixture (`3comp`): butan-1-ol / propan-2-ol / water
- five component mixture (`5comp`): methanol / ethanol / propan-2-ol / butan-1-ol / acetone

Each mixture is run with five liquid-phase models (NRTL, modified UNIFAC, COGNAC explicit, HANNA without ensemble, HANNA). Pure-component properties come from [`thermo_database/`](thermo_database). The column starts at total reflux and runs until `tmax` or reboiler dry-out. The script prints the solve time per model and stores the reboiler and reflux-drum compositions in `results`.

## Running it

```julia
using Pkg
Pkg.activate("examples/cognac_comparison")
Pkg.instantiate()
include("examples/cognac_comparison/run_comparison.jl")
```
