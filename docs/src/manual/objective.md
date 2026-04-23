# [Objective function](@id man_objective)
The basic objective for computing a control action in the MPC controller is of the form
```math
\sum_{k=0}^{N-1} {\left((Cx_{k}-r)^T Q (C x_{k}-r) + u_{k}^T R u_{k} + \Delta u_{k}^T R_r \Delta u_k + l_k^T u_k\right)},
```
where $N$ is the prediction horizon of the controller.

For an MPC controller `mpc`, the weighting matrices $Q$, $R$, and $R_r$ that defines the objective can be set with 

```@tab
# julia
set_objective!(mpc;Q,R,Rr)
# python
mpc.set_objective(Q=Q, R=R, Rr=Rr)
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

```@tab
# julia
mpc.settings.reference_preview = true
setup!(mpc)
# python
mpc.settings({"reference_preview": True})
mpc.setup()
```

Then provide a reference trajectory matrix of size `(ny, Np)` to `compute_control`:

```@tab
# julia
r_trajectory = [1.0 1.5 2.0 2.0 2.0;   # Reference for output 1
                0.0 0.0 0.5 1.0 1.0]   # Reference for output 2
u = compute_control(mpc, x; r=r_trajectory)
# python
import numpy as np
r_trajectory = np.array([[1.0, 1.5, 2.0, 2.0, 2.0],   # Reference for output 1
                          [0.0, 0.0, 0.5, 1.0, 1.0]])  # Reference for output 2
u = mpc.compute_control(x, r=r_trajectory)
```

## Linear control cost
The term $l_k^T u_k$ adds an optional linear cost on the control signal. This is useful in economic MPC where control actions have time-varying costs (e.g., electricity prices).

Enable linear cost with:

```@tab
# julia
mpc.settings.linear_cost = true
setup!(mpc)
# python
mpc.settings({"linear_cost": True})
mpc.setup()
```

Then provide the linear cost to `compute_control`:

```@tab
# julia
# Constant linear cost across horizon
u = compute_control(mpc, x; l=[0.5])

# Time-varying linear cost (nu × Np matrix)
l_trajectory = [0.1 0.2 0.5 0.8 1.0]  # Cost varies over prediction horizon
u = compute_control(mpc, x; l=l_trajectory)
# python
import numpy as np
# Constant linear cost across horizon
u = mpc.compute_control(x, l=[0.5])

# Time-varying linear cost (nu x Np array)
l_trajectory = np.array([[0.1, 0.2, 0.5, 0.8, 1.0]])  # Cost varies over prediction horizon
u = mpc.compute_control(x, l=l_trajectory)
```

## Affine parameters in objective and constraints
LinearMPC also supports a stagewise affine parameter trajectory $p_k$ entering the problem as
```math
(E p_k + e)^T u_k
```
in the objective, and through additional constraint terms such as
```math
\underline{b} \le A_x x_k + A_u u_k + A_p p_k \le \overline{b}.
```

The matrices `E` and `Ap` are stagewise coefficients. A constant parameter vector is broadcast across the prediction horizon, while an `(np, Np)` matrix gives a time-varying parameter preview.

```@tab
# julia
set_objective!(mpc; Q=[1.0], R=[0.1], E=[1.0 0.0], e=[0.2])
add_constraint!(mpc; Au=[1.0;;], Ap=[0.5 1.0], ub=[1.0], ks=1:mpc.Np)

u = compute_control(mpc, x; p=[0.3, -0.1])

# Time-varying parameter preview
p_preview = [0.3 0.2 0.1 0.1;
             -0.1 -0.1 0.0 0.1]
u = compute_control(mpc, x; p=p_preview)
# python
mpc.set_objective(Q=[1.0], R=[0.1], E=[[1.0, 0.0]], e=[0.2])
mpc.add_constraint(Au=[[1.0]], Ap=[[0.5, 1.0]], ub=[1.0], ks=range(1, mpc.Np + 1))

u = mpc.compute_control(x, p=[0.3, -0.1])

# Time-varying parameter preview
p_preview = np.array([[0.3, 0.2, 0.1, 0.1],
                      [-0.1, -0.1, 0.0, 0.1]])
u = mpc.compute_control(x, p=p_preview)
```
