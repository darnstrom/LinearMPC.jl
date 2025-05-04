# [Models](@id man_model)
As its name suggests **Model** Predictive Control uses a predictive model when computing a control action. The models used internally in **LinearMPC.jl** are **discrete-time state-space models** of the form 
```math
x_{k+1} = F x_k + G u_k
```
where the state at time step $k+1$ is a linear combinations of the current state $x_k$ and a control action $u_k$. If matrices $F$ and $G$ that define such a state-space model is available, a MPC controller that uses that model is created with 

```julia
mpc = LinearMPC.MPC(F,G)
```
One can also create a `Model` struct first, and then create the MPC controller, with 

```julia
model = LinearMPC.Model(F,G)
mpc = LinearMPC.MPC(model)
```

## Continuous-time models 
Often, the system dynamics is given in continuous time rather than discrete time. That is, the system dynamics is a state space model of the form 

```math
\frac{d}{dt} x(t) = A x(t) + B u(t).

```
One alternative to deal with this is to [discretize](https://en.wikipedia.org/wiki/Discretization#discrete_function) the system, to yield matrices $F$ and $G$ instead of $A$ and $B$. A MPC controller that use a discretized continuous-time model can be created with

```julia
mpc = LinearMPC.MPC(A,B,Ts)
```
where `Ts` is the _sample time_, which  and determines how much one time step is. (For example, $T_s = 0.1$ corresponds to 10 time steps being 1 second.) As for discrete-time state spaces models, one can first create a `Model` struct, and create a MPC controller based on this model as follows: 

```julia
model = LinearMPC.Model(A,B,Ts)
mpc = LinearMPC.MPC(model)
```


## Outputs
In the MPC, one can make an "output" $y$ of the system track a reference $r$. This "output" is a linear combination of the form  
```math
y_k = C x_k, 
```
where $C$ is a matrix that maps the states output. 

!!! note "Difference to measurements"
    In this context, $y$ is not necessarily what is measured from the system. Instead, it should be thought of quantities that are of interest to be controlled.

The output of the system can be set with the  optional argument `C`. So for the model 
```math
x_{k+1} = F x_k + G u_k, \quad y_k = C x_k,
```
an MPC controller that uses this model can be created with

```julia
mpc = LinearMPC.MPC(F,G;C)
```
or by first creating a `Model` as 
```julia
model = LinearMPC.Model(F,G;C)
mpc = LinearMPC.MPC(model)
```

Similarly, an MPC controller of the continuous-time system
```math
\frac{d}{dt} x(t) = A x(t) + B u(t), \quad y(t) = C x(t)
```
is setup with 
```julia
mpc = LinearMPC.MPC(A,B,Ts;C)
```
or 
```julia
model = LinearMPC.Model(A,B,Ts;C)
mpc = LinearMPC.MPC(model)
```

!!! note "Default value"
    If $C$ is not provided, it is assumed that all of the states are the output (that is, $y=x$, or, equivalently, $C=I$.)
## Disturbances
Models with disturbances $d$ are also supported. Specifically, linear combinations of $d$ is allowed to be additive to the dynamics and to the measurements as  
```math
x_{k+1} = F x_k + G u_k + G_d d_k, \quad y_k = C x_k + D_d d_k, 
```
where the matrices $G_d$ and $D_d$ determines the linear combination for the dynamics and the output, respectively. Both $G_d$ and $D_d$ can be set as optional arguments:

```julia
mpc = LinearMPC.MPC(F,G;Gd,C,Dd)
```
Similarly, continuous-time state space models of the form 
```math
\frac{d}{dt} x(t) = A x(t) + B u(t) + B_d d(t), \quad y(t) = C x(t) + D_d x(t)
```
can be used with 
```julia
mpc = LinearMPC.MPC(A,B,Ts;Bd,C,Dd)
```

## Linearization
For a nonlinear discrete-time state-space model of the form  
```math
x_{k+1} = f(x_k,u_k,d_k), \quad y_k = h(x_k,u_k,d_k),    
```
an MPC controller with a linearized model around an operating point $x_o$, $u_o$ and $d_o$ can be created with 
```julia
model = LinearMPC.Model(f,h,x_o,u_o; d=d_o)
mpc = LinearMPC.MPC(model)
```

Similary, for a nonlinear _continuous-time_ state-space model of the form

```math
\frac{d}{dt} x(t) = f(x(t),u(t),d(t)), \quad y(t) = h(x(t),u(t),d(t))
```
an MPC controller with a linearized model around an operating point $x_o$, $u_o$ and $d_o$ can be created with 
```julia
model = LinearMPC.Model(f,h,x_o,u_o,Ts; d=d_o)
mpc = LinearMPC.MPC(model)
```

## Using ControlSystems.jl
**Linear MPC** also support system models from [ControlSystems.jl](https://juliacontrol.github.io/ControlSystems.jl/stable/man/creating_systems/). Specifically it supports these types of systems: 

* Transfer functions
* State-space models
* Delayed LTI systems

For a system `sys` created with ControlSystems.jl, an MPC controller can simply created with  
```julia
mpc = LinearMPC.MPC(sys)
```
For example, given the transfer function  
```math
 G(s) = \frac{1}{s^2+2s+1} 
```
and MPC controller using this transfer function as its model can be created with  
```julia
sys = tf([1.0],[1,2,1])
mpc = LinearMPC.MPC(sys)
```
