# [Simulating](@id man_simulation)

`Simulation` runs a closed-loop simulation of a controller and stores the signals in a
plot-friendly object.

```julia
sim = Simulation(mpc; x0 = x0, N = 100, r = r)
```

The returned object contains

- `sim.ts`: time vector
- `sim.xs`: simulated states
- `sim.ys`: controller outputs
- `sim.us`: applied controls
- `sim.rs`: references used for plotting
- `sim.ds`: disturbances used by the simulated plant
- `sim.xhats`: observer state estimates
- `sim.yms`: measured outputs passed to the observer
- `sim.solve_times`: time spent in `compute_control`

If a plotting backend is loaded, the built-in recipe can be used directly:

```julia
using Plots
plot(sim)
```

## Custom Dynamics

By default, `Simulation` uses the model's `true_dynamics`. A custom plant can be passed
as the first argument:

```julia
dynamics = (x, u, d) -> Aplant * x + Bplant * u + Ed * d
sim = Simulation(dynamics, mpc; x0 = x0, N = 100, r = r)
```

The measurement function can also be overridden:

```julia
get_measurement = (x, d) -> Cmeas * x + noise()
sim = Simulation(dynamics, mpc; x0 = x0, get_measurement)
```

## Reference, Disturbance, and Parameter Trajectories

`r`, `d`, and `p` may be vectors or matrices. If a trajectory matrix is shorter than
the simulation horizon, the final column is held.

When preview settings are enabled, `Simulation` automatically extracts horizon previews
before each control solve:

```julia
mpc.settings.reference_preview = true
mpc.settings.disturbance_preview = true
mpc.settings.parameter_preview = true

sim = Simulation(mpc; x0 = x0, N = 100, r = r_traj, d = d_traj, p = p_traj)
```

## Dynamic Preview Callback

Some controllers need references, disturbances, or generalized parameters that depend on
the current observer estimate. For this, `Simulation` accepts a `preview` callback. The
callback is evaluated after observer correction and before `compute_control`:

```julia
preview = function (mpc, xhat, y, k)
    r_preview = make_reference_preview(xhat, k)
    p_preview = make_parameter_preview(xhat, k)
    return (r = r_preview, p = p_preview)
end

sim = Simulation(dynamics, mpc; x0 = x0, N = 100, preview)
```

The callback can return any subset of `r`, `d`, and `p`. Values returned by the callback
override the static trajectories for the current step. If `r` is returned, the first
column is also stored in `sim.rs` so that automatic plotting shows the dynamic reference.

This is useful for offset-free periodic tracking. A typical pattern is:

```julia
mpc.settings.reference_preview = true
mpc.settings.disturbance_preview = true
mpc.settings.parameter_preview = true
set_objective!(mpc; Q, Qf, R, Eu = -R)

set_offset_free_observer!(mpc;
    method = :periodic,
    Bd = Bd,
    Cd = Cd,
    period = period,
    Q = Qobs,
    R = Robs,
)

preview = function (mpc, xhat, y, k)
    xbar, ubar = periodic_offset_free_target_preview(
        mpc.state_observer, F, G, Cz, reference_period(k);
        Bd = Bd,
        Cd = Cd,
        Np = mpc.Np,
    )
    return (r = xbar[:, 2:end], p = ubar)
end

sim = Simulation(nonlinear_dynamics, mpc; x0 = x0, N = 1000, preview)
```

Here `xbar[:, 2:end]` is the state-reference preview for the predicted states and
`ubar` enters through the generalized-parameter input cost. The disturbance preview is
assembled automatically from the periodic offset-free observer.

## Callback

A separate `callback` can be used for logging or side effects after each simulated plant
step:

```julia
callback = (x, u, d, k) -> println("step ", k, ": ", x)
sim = Simulation(mpc; x0 = x0, N = 100, callback)
```
