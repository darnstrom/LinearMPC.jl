# [Objective function](@id man_objective)
The basic objective for computing a control action in the MPC controller is of the form 
```math
\sum_{k=0}^{N-1} {\left((Cx_{k}-r)^T Q (C x_{k}-r) + u_{k}^T R u_{k} + \Delta u_{k}^T R_r \Delta u_k\right)},
```
where $N$ is the prediction horizon of the controller.

For an MPC controller `mpc`, the weighting matrices $Q$, $R$, and $R_r$ that defines the objective can be set with 
```julia
set_objective!(mpc;Q,R,Rr)
```
If a vector value of $Q$,$R$, or $R_r$ are inputted to `set_objective`, it will be interpreted as a diagonal matrix with the vector on diagonal. Similarly, if a scalar is inputted, it will be interpreted as a diagonal matrix with this value on the diagonal.

The default values are $Q=I$, $R=I$, and $R_r = 0$.

## Reference tracking
The first term $(Cx_{k}-r)^T Q (C x_{k}-r)$ is for _reference tracking_, which means that it aims to make $Cx = r$, where the matrix $C$ is set by the user when initializing the model. Typically, the weighting matrix $Q$ is a diagonal matrix, where higher weights are selected for the rows of $C x$ that should be prioritized.

## Control energy
The second term $u_k^T R u_k$ penalizes larger values of the control action. Typically, the weighting matrix $R$ is a diagonal matrix, where higher weights are selected for the control signals that are more "expensive".

## Control change 
The third term $\Delta u_k^T R \Delta u_k$ penalizes changes to the control action. Typically, the weighting matrix $R_r$ is a diagonal matrix, where higher weights are selected for the control signals that are more "expensive".

## Terminal constraint
A term of the form $(Cx_N-r)^T Q_f (C x_N -r)^T$ can also be added to the objective. This adds an additional terminal cost, which can be useful to ensure stability[^Mayne00]. The terminal weight $Q_f$ is set with `set_objective!`, and is by default the same as $Q$. 


[^Mayne00]: Mayne, David Q., et al. "Constrained model predictive control: Stability and optimality." _Automatica_ 36.6 (2000): 789-814.

## Cross term 
It is also possible to include a cross term $x_k^T S u_k$ in the objective. This term can also be set with the `set_objective!` function. By default, $S=0$.

## Reference Preview
By default, the reference tracking term uses a constant reference $r$ across the prediction horizon: $(Cx_{k}-r)^T Q (C x_{k}-r)$ for all $k$. 

Reference preview allows time-varying references $r_k$ in the objective:
```math
\sum_{k=0}^{N-1} (Cx_{k}-r_k)^T Q (C x_{k}-r_k)
```

Enable reference preview with:
```julia
mpc.settings.reference_preview = true
setup!(mpc)
```

Then provide a reference trajectory matrix of size `(ny, Np)` to `compute_control`:
```julia
r_trajectory = [1.0 1.5 2.0 2.0 2.0;   # Reference for output 1
                0.0 0.0 0.5 1.0 1.0]   # Reference for output 2  
u = compute_control(mpc, x; r=r_trajectory)
```

