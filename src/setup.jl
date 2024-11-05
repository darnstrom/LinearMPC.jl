function update_dynamics!(mpc::MPC,F,G)
    mpc.F,mpc.G = F,G
    mpc.mpQP = mpc2mpqp(mpc)
    return
end
