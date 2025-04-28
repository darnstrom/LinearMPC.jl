# [Example - Inverted Pendulum on a Cart](@id ex_invpned)
To exemplify the basics of **LinearMPC.jl**, we consider the classical example of controlling a inverted pendulum on a cart. The cart can be moved by applying a force $u$, to control the cart position $p$ 
and pendulum angle $\theta$. The start of using MPC is to have a _model_ of the system we want to control, and the model which is used in an MPC controller is discrete-time state-space models of the form $x_{k+1} = F x_k G u_k$. In ref we discuss more how to use different types of models, but the end-point is always a discrete-time state-space system. For simplicity, we assume that the following state-space model is given for cart system 

```m̀ath
F = ,\quad G = 
```
where $p$ and $\dot{p}$ are the position/velocity of the cart, and $\theta$ and $\dot{theta}$ are the angle and angular velocity of the pendulum. 

Our goal is to control the cart's position $p$ and the angle of the pendulum $\theta$. We denote the signals we want to control with $y$, which in our case gives $y_1 = p$ and $y_2 = \theta$. How the states relates to our output $y$ can compactly be written as $y = Cx$, which for the cart system gives 
```math
C = 
```
An additional requirement is that we don't want the control signal to be too "jittery". This is addressed with a penalty $\Delta u R_r \Delta u$, where $\Delta u$ denotes the change in the control input between to time steps. 
The objective at a single time step, thus becomes $(Cx - y)^T Q (Cx-y) + \Delta u^T R \Delta u$.


In addition to the dynamics of the system, we also want to impose additional constraint when controlling the system. First, there is a limited amount of force $u$ that can be applied. Specifically, this gives us the following bounds on the control: $-2 \leq u \leq 2$. Moreover, we also want the angle $\theta$ to not be to large, since this would our model. Therefore, we also have the constraint $-0.2 \leq \theta \leq 0.2$.




```julia
F = 
```

Taken together, we pose the control problem as a receding horizon control problem. 

Control horizon / Prediction horizon (hint at move blocking)

Simulate the system

## Generating C-code 
