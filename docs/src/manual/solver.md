# [Solver Settings](@id man_solver)

LinearMPC uses the [DAQP](https://darnstrom.github.io/daqp/) solver for solving quadratic programs. You can configure solver settings to adjust tolerances, iteration limits, and other parameters.

## Changing Settings

After setting up the MPC controller, you can modify solver settings using `DAQP.settings`:

```julia
using LinearMPC

# Create and setup MPC
mpc = MPC(A, B; Np=10)
set_bounds!(mpc; umin=-1, umax=1)
setup!(mpc)

# Change solver settings
DAQP.settings(mpc.opt_model, Dict(
    :iter_limit => 2000,
    :primal_tol => 1e-8
))

# Compute control with new settings
u = compute_control(mpc, x)
```

## Basic Settings

For full documentation of all DAQP settings, see the [DAQP settings reference](https://darnstrom.github.io/daqp/parameters/#settings).


|  Parameter |  Description| Default |
|:-------------|:------------------|:------:|
| `primal_tol`  | Tolerance for primal infeasibility|  1e-6 |
| `dual_tol`	 | Tolerance for dual infeasibility| 1e-12|
| `progress_tol` | Minimum change in objective function to consider it progress | 1e-6|
| `cycle_tol` | Allowed number of iterations without progress before terminating| 10 |
| `iter_limit` | Maximum number of iterations before terminating| 1000 |
| `rho_soft` | Weight used for soft constraints (higher enables more violations) | 1e-6|

## Code Generation

When generating C code, you can pass solver settings via the `opt_settings` argument:

```julia
codegen(mpc; opt_settings=Dict(:iter_limit => 500))
```

These settings will be embedded in the generated C code.

