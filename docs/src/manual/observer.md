# [State observer](@id man_observer)

Model predictive control requires _all_ the state of the system to be known. All states are, however, rarely measured in practice. Therefore, a _state observer_ that estimates the current states based on measurements is needed. The goal of a state observer is therefore to generate an estimate $\hat{x}$ of the true state $x$, based on measurements $y$ and applied controls $u$ such that $\hat{x} \approx x$

**LinearMPC.jl** provides a simple steady-state Kalman filter for estimating the states. The steady-state Kalman filter assumes the following model of the dynamical system
```math
x_{k+1} = F x_k + G u_k + w_k, \qquad y_k = C x_k + e_k
```
where the process noise $w_k \sim \mathcal{N}(0,Q)$ and the measurement noise $e_k \sim \mathcal{N}(0,R)$. Each time instance, the filter performs two steps: a _correction step_ and a _prediction step_. 
In the correction steps, the current state estimate $\hat{x}$ is corrected based on a measurement $y$ according to   
```math
\hat{x} \leftarrow \hat{x} + K(y-C \hat{x}),
```
which corrects the state based on the deviation between the measurement $y$ and the _expected_ measurement $C\hat{x}$.
The correction step is followed by a prediction step, where the state at the next time step is predicted based on the applied control $u$:
```math
\hat{x} \leftarrow F \hat{x} + G u.
```

For a MPC struct `mpc`, a steady state Kalman filter can be created with
```julia  
set_state_observer!(mpc;Q,R)
```
where `Q` and `R` are covariance matrices for the process noise and measurement noise, respectively. These are used as tuning variables, where the relative size of the elements in `Q` and `R` determines if the prediction or the measurements should be trusted more. By default, `set_state_observer!` uses $F$, $G$, and $C$ from the MPC structure. It is possible to override this by passing the optional arguments `F`,`G`,`C` to `set_state_observer!`.

## Getting, setting, correcting, and predicting state
The current state of the observer is accessed with 
```julia
x = get_state(mpc)
```

The state of the observer can be set to `x0` with
```julia
set_state!(mpc,x0)
```

Given a measurement `y`, the state of the observer can be corrected with
```julia
correct_state!(mpc,y)
```

Given a control action `u`, the next state can be predicted with 
```julia
predict_state!(mpc,u)
```

## Using the observer  
Typically, you want to do a correction step before computing a new control action, and then use the new control action to predict the next state. A typical sequence of calls are therefore 
```julia
x = correct_state!(mpc,y)
u = compute_control(mpc,x)
predict_state!(mpc,u)
```

## Code generation
If an observer has been set, the `codegen` function will in addition to C-code for the MPC also generate C-code for the observer. Specifically the C functions `mpc_correct_state(state, measurement, disturbance)` and  `mpc_predict_state(state, control, disturbance)` will be generated, where `state`, `measurement`, `control` and `disturbance` are floating-point arrays. The updated state will be available in the floating-point array `state` after calling the functions. 

Similar to the previous section, the following is a typical sequence of calls when using an observer together with `mpc_compute_control` 

```c
mpc_correct_state(state,measurement,disturbance)
mpc_compute_control(control,state,reference,disturbance)
mpc_predict_state(state,control,disturbance)
```
Note for cases when there is no distrubance, one can pass `NULL` instead of `disturbance` as last argument to the functions.
