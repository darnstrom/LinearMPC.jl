# Benchmark
Code for producing benchmark results is available [here](https://github.com/darnstrom/lmpc-codegen-benchmark).

Here show a closed-loop benchmark for a inverted-pendulum on a cart example, sweeping the prediction and control horizons together over `N = 50, 75, 100, 125`. Both input and state constraints are imposed. 

Specifically, similar results can be genereated with the following commands.

```bash
uv sync
uv run python run_all.py \
  --problem inverted_pendulum \
  --scaling \
  --solvers lmpc casadi cvxpygen acados tinympc \
  --horizons 50 75 100 125 \
  --results-dir results
```


In this benchmark, **LinearMPC.jl** (denoted `lmpc` in the plots since its sister Python package [lmpc](https://github.com/darnstrom/lmpc) is used) produces the fastest code across the considered horizons. Moreover, **LinearMPC.jl** has the smallest memory footprint. Note, however, that the dense formulation that **LinearMPC.jl** will sooner or later exceed that of the sparse solvers, as can been seen by the slope memory footprint.

```@raw html
<div style="text-align: center;">
    <img src="../../assets/benchmark_scaling_time.png" style="width: 70%; height: auto;">
    <br>
    <img src="../../assets/benchmark_scaling_memory.png" style="width: 70%; height: auto;">
</div>
```


**Detailed information for horizon 50:**

```@raw html
<div style="text-align: center;">
    <img src="../../assets/benchmark_timing_cdf.png" style="width: 70%; height: auto;">
    <br>
    <img src="../../assets/benchmark_memory.png" style="width: 70%; height: auto;">
</div>
```

Note that the memory footprint for `acados` is large due to its dependence on `blasfeo` and `hpipm`. The memory footprint could be reduced if a more targeted compilation is performed.
