# [Code Generation](@id man_codegen)

Most real-time controllers run on embedded hardware, which often require the controller to be implemented in a low-level programming language like C. However, implementing an MPC controller in C from scratch is a very time consuming endeavor. To simplify the process, **LinearMPC.jl** can generate C-code for MPC controllers that have been designed and tested in Julia, which enables the MPC controller to easily be applied on embedded systems. To generate such C code in a directory `code_dir` for an MPC controller `mpc`, we can run the following code

```@tab
# julia
LinearMPC.codegen(mpc; dir="code_dir", fname="test_mpc")
# python
mpc.codegen(fname="test_mpc", dir="code_dir")
```

where `fname` determines some of the naming of the generated code.

The main function of interest (located in `{fname}.h`) is `mpc_compute_control(control, state, reference, disturbance)`. This function computes the optimal control given the current `state`, `reference`, and measured disturbances `disturbance`, which are all floating-point arrays. The optimal control is stored in the floating-point array `control`.

!!! note "Previous control action"
    If there is a penalty on the change in of control actions $\Delta u$, the function `mpc_compute_control` uses the value that is in `control` as the previous control action `uprev`. 

How the values of `state`,`reference`, and `disturbance` in the generated C-code are set depends on the particular application. For example, `state` might come from a state observer, and `reference` might come from some motion planner or user interface.
