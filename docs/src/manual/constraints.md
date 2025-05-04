# [Constraints](@id man_constraints)
The constraints that are supported in **LinearMPC.jl** are of the form 
```math
\underline{b} \leq A_x x_k + A_u u_k \leq \overline{b}.
```
Such constraints can be added with the function `add_constraint!`.

Bounds on controls ($\underline{u} \leq u \leq \overline{u}$) and on outputs ($\underline{y} \leq y \leq \overline{y}$) can be set in a special way with the function `set_bounds!`.


## Illustrative example
As an illustrative example, consider the case with two controls $u_1$ and $u_2$, two outputs $y_1$ and $y_2$, and three states $x_1$, $x_2$, and $x_3$, where we want to impose the following constraints
```math
\begin{aligned}
-2 &\leq u_1 &\leq 1  \\
0 &\leq u_2 &\leq 3 \\
-1 &\leq y_1 &\leq 1  \\
-3 &\leq y_2 &\leq 0 \\
-1 &\leq x_2 - 3 x_3 &\leq 5 \\
-1 &\leq x_1 - 2 x_2 + u_2 &\leq 1
\end{aligned}
```
This can be added to an MPC controller `mpc` with the following commands
```julia
# control + output bounds
set_bounds!(mpc;umin=[-2,0], umax = [1,3], ymin = [-1,-3], ymax = [1,0])

# -1 ≤ x2 - 3 x3 ≤  5 
add_constraint!(mpc; Ax = [1 0 3], lb = -2, ub = 5)

# -1 ≤ x1 - 2 x2 + u2 ≤  1
add_constraint!(mpc; Ax = [1 -2 0], Au = [0 1], lb = -1, ub = 1)
```

## Constraint horizon 
The constraints from above will be enforced for all time steps up until the prediction horizon. It is, however, possible to specify a limited subset of time steps when a constraint should be satisfied. In other words, constraints are of the form
```math
\underline{b} \leq A_x x_k + A_u u_k \leq \overline{b},\qquad  \text{for all } k\in \mathcal{K},
```
where the index set ${\mathcal{K}} \subseteq \{1,\dots,N_p\}$ determines for which time steps the constraint should be enforced. By default $\mathcal{K} = \{1,\dots N_p\}$. The indices for which the constraints should be enforced are specified with the optional argument `ks` to `add_constraint!`.  

For example, say that we want to add the constraint

```math
-1 \leq x_1 + x_2 + 3 x_3 \leq 1, \qquad k \in \{2,4,6\}
```

to the MPC controller `mpc`. This can be done with 
```julia
add_constraint!(mpc;Ax=[1 1 3], lb = -1, ub = 1, ks = [2,4,6])
```



## Soft constraints
To ensure that there is a solution to the resulting optimization problem, it can be good to _soften_ some constraint. A soft constraints means that the constraints is allowed to be violated, but such violations are penalized. A constraint can be marked to be soft by setting the optional argument `soft` to `true`. For example, if the constraint from above should be soft, one could run the command 
```julia
add_constraint!(mpc;Ax=[1 1 3], lb = -1, ub = 1, ks = [2,4,6], soft=true)
```
By default, bound constraints on inputs $u$ are _hard_, and bound constraint outputs $y$ are _soft_.

To be more specific, a constraint $\underline{b} \leq A_x x_k + A_u u_k \leq \overline{b}$ is softened by introducing a slack variable $\epsilon$, and replace the constraint with the constraints 
```math
\begin{aligned}
\underline{b}-\epsilon \leq &A_x x_k + A_u u_k & \\
& A_x x_k + A_u u_k &\leq \overline{b} + \epsilon.
\end{aligned}
```
To penalize violation of the constraint, a term $\rho \|\epsilon\|$ is added to the objective function, where the weight $\rho$ determines how hard the soft constraints should be penalized. By default, $\rho$ is set to `1e6`. 
