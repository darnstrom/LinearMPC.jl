struct MPCExample
    name::String
    mpc::MPC
    range::ParameterRange
    scenarios::Vector{Scenario}
end

MPCExample(name::AbstractString, mpc::MPC, range::ParameterRange; scenarios=Scenario[]) =
    MPCExample(String(name), mpc, range, Scenario[scenarios...])

_normalize_example_id(s::Union{AbstractString,Symbol}) = replace(lowercase(String(s)), r"[\s_-]" => "")

const _MPC_EXAMPLE_SPECS = (
    invpend = (
        name = "Inverted Pendulum on a Cart",
        aliases = ("inv_pend", "invpend", "invpendcart"),
        defaults = (Np = 50, Nc = 5),
    ),
    dcmotor = (
        name = "DC Motor Position Control",
        aliases = ("dc_motor", "dcmotor"),
        defaults = (Np = 10, Nc = 2),
    ),
    aircraft = (
        name = "Aircraft",
        aliases = ("aircraft",),
        defaults = (Np = 10, Nc = 2),
    ),
    chained = (
        name = "Chained",
        aliases = ("chained", "chained-firstorder"),
        defaults = (Np = 10, Nc = 10),
    ),
    mass_spring = (
        name = "Mass-Spring System",
        aliases = ("mass-spring", "mass", "spring"),
        defaults = (Np = 10, Nc = 10),
    ),
    nonlinear_demo = (
        name = "Linearized Nonlinear Demo",
        aliases = ("nonlinear", "nonlin"),
        defaults = (Np = 5, Nc = 2),
    ),
    invpend_contact = (
        name = "Colliding Inverted Pendulum on a Cart",
        aliases = ("invpend_contact",),
        defaults = (Np = 10, Nc = 10),
    ),
    ballplate = (
        name = "Ball and Plate",
        aliases = ("ball", "ballplate"),
        defaults = (Np = 10, Nc = 2),
    ),
    crazyflie = (
        name = "Quadcopter",
        aliases = ("quad", "quadcopter", "crazyflie"),
        defaults = (Np = 10, Nc = 10),
    ),
    satellite = (
        name = "Satellite",
        aliases = ("satellite",),
        defaults = (Np = 20, Nc = 20),
    ),
    rocket = (
        name = "Rocket",
        aliases = ("rocket",),
        defaults = (Np = 100, Nc = 10),
    ),
)

const _MPC_EXAMPLE_ALIASES = let aliases = Dict{String,Symbol}()
    for (id, spec) in pairs(_MPC_EXAMPLE_SPECS)
        aliases[_normalize_example_id(String(id))] = id
        for alias in spec.aliases
            aliases[_normalize_example_id(alias)] = id
        end
    end
    aliases
end

_mpc_example_ids() = collect(propertynames(_MPC_EXAMPLE_SPECS))
mpc_example_names() = sort!(String.(collect(propertynames(_MPC_EXAMPLE_SPECS))))

function _resolve_example_id(s::Union{AbstractString,Symbol})
    id = get(_MPC_EXAMPLE_ALIASES, _normalize_example_id(s), nothing)
    isnothing(id) || return id
    error("There is currently no example `$(s)`. Available examples: $(join(mpc_example_names(), ", ")).")
end

_default_horizons(id::Symbol) = begin
    defaults = getfield(_MPC_EXAMPLE_SPECS, id).defaults
    return defaults.Np, defaults.Nc
end

_finalize_example(id::Symbol, mpc::MPC, range::ParameterRange; scenarios=Scenario[]) =
    MPCExample(getfield(_MPC_EXAMPLE_SPECS, id).name, mpc, range; scenarios)

function _merge_example_kwargs(params, kwargs)
    merged = Dict{Symbol,Any}()
    for (key, value) in pairs(params)
        merged[Symbol(key)] = value
    end
    for (key, value) in pairs(kwargs)
        merged[key] = value
    end
    return merged
end

function _build_mpc_example(::Val{:invpend}, Np, Nc; settings=nothing, kwargs...)
    M = get(kwargs, :M, 1.0)
    m = get(kwargs, :m, 1.0)
    l = get(kwargs, :l, 0.5)
    damp = get(kwargs, :damp, 10.0)

    g = 9.81
    scale, Mm = 100, M + m

    f = (x, u, d) -> [x[2];
                      (scale*u[1] - damp*x[2] - m*l*x[4]^2*sin(x[3]) + m*g*sin(x[3])*cos(x[3])) / (M + m*sin(x[3])^2);
                      x[4];
                      (g*sin(x[3]) + (scale*u[1] - damp*x[2] - m*l*x[4]^2*sin(x[3]))*cos(x[3]) / Mm) / (l - m*l*cos(x[3])^2 / Mm)]
    h = (x, u, d) -> [x[1]; x[3]]

    Ts = 0.01
    xo, uo = zeros(4), zeros(1)
    model = LinearMPC.Model(f, h, xo, uo, Ts)

    mpc = MPC(model; Nc, Np)
    set_objective!(mpc; Q = [1.2^2, 1], R = [0.0], Rr = [1.0])
    set_bounds!(mpc; umin = [-2.0], umax = [2.0])
    isnothing(settings) || (mpc.settings = settings)

    range = ParameterRange(mpc)
    range.xmax[:] .= 20 * ones(4)
    range.xmin[:] .= -20 * ones(4)
    range.rmax[:] .= 20 * ones(2)
    range.rmin[:] .= -20 * ones(2)
    range.dmax[:] .= 20 * ones(1)
    range.dmin[:] .= -20 * ones(1)

    scenarios = [
        Scenario([0.0, 0.0, 0.15, 0.0]; T = 2.0, r = [0.0, 0.0]),
        Scenario(zeros(4); T = 2.0, r = [1.0, 0.0]),
    ]
    return _finalize_example(:invpend, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:dcmotor}, Np, Nc; settings=nothing, kwargs...)
    A = [0 1.0 0 0; -51.21 -1 2.56 0; 0 0 0 1; 128 0 -6.401 -10.2]
    B = 440 * [0; 0; 0; 1.0;;]
    C = [1 0 0 0; 1280 0 -64.01 0]
    Ts = 0.1
    tau = 78.5398
    C = C ./ [2*pi; 2*tau]

    mpc = MPC(A, B, Ts; C, Np, Nc)
    set_objective!(mpc; Q = [0.1^2, 0], R = [0.0], Rr = [0.1^2])
    set_bounds!(mpc; umin = [-0.5], umax = [0.5])
    add_constraint!(mpc, Ax = C[2:2, :], lb = [-0.5], ub = [0.5], ks = 2:min(mpc.Nc + 2, mpc.Np), soft = true)

    if isnothing(settings)
        mpc.settings.reference_tracking = true
    else
        mpc.settings = settings
    end

    range = ParameterRange(mpc)
    range.xmax[:] = [4*pi 4*pi/Ts 4*pi*20 4*pi*20/Ts]
    range.xmin[:] = -[4*pi 4*pi/Ts 4*pi*20 4*pi*20/Ts]
    range.rmax[:] .= [5; 0.5]
    range.rmin[:] .= -[5; 0.5]
    range.umax[:] .= [0.5023]
    range.umin[:] .= -[0.5023]

    scenarios = [Scenario(zeros(4); T = 2.0, r = [1.0, 0.0])]
    return _finalize_example(:dcmotor, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:aircraft}, Np, Nc; settings=nothing, kwargs...)
    A = [-0.0151 -60.5651 0 -32.174;
         -0.0001 -1.3411 0.9929 0;
         0.00018 43.2541 -0.86939 0;
         0 0 1 0]
    B = [-2.516 -13.136; -0.1689 -0.2514; -17.251 -1.5766; 0 0]
    C = [0 1.0 0 0; 0 0 0 1]

    Ts = 0.05
    F, G = zoh(A, B, Ts)
    C = C ./ [1; 200]
    Dd = [1.0 0; 0 200] ./ [1; 200]

    mpc = MPC(F, 50 * G; C, Np, Nc, Ts, Dd)
    set_objective!(mpc; Q = [10, 10].^2, R = zeros(2), Rr = [0.1, 0.1].^2)
    set_bounds!(mpc, umin = [-0.5; -0.5], umax = [0.5; 0.5])
    set_output_bounds!(mpc, ymin = [-0.5; -0.5], ymax = [0.5; 0.5], ks = 2:2)

    if isnothing(settings)
        mpc.settings.reference_tracking = true
    else
        mpc.settings = settings
    end

    range = ParameterRange(mpc)
    range.xmax[:] .= 20 * ones(4)
    range.xmin[:] .= -20 * ones(4)
    range.dmax[:] .= 20 * ones(2)
    range.dmin[:] .= -20 * ones(2)
    range.rmax[:] .= [1; 0.05]
    range.rmin[:] .= -[1; 0.05]

    scenarios = [Scenario(zeros(4); T = 2.0, r = [0.1, 0.0], d = zeros(2))]
    return _finalize_example(:aircraft, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:chained}, Np, Nc; settings=nothing, kwargs...)
    nx = Int(get(kwargs, :nx, 1))
    A = -Matrix(I, nx, nx) + diagm(-1 => ones(nx - 1))
    B = [1; zeros(nx - 1, 1);;]
    C = Matrix(I, nx, nx)
    Ts = 1
    F, G = zoh(A, B, Ts)

    mpc = MPC(F, G; C, Np, Nc, Ts)
    set_objective!(mpc; Q = ones(nx), R = [0.0], Rr = [1.0])
    set_bounds!(mpc, umin = [-1.0], umax = [1.0])
    set_output_bounds!(mpc, ymin = -10 * ones(nx), ymax = 10 * ones(nx), ks = 2:mpc.Nc)

    if isnothing(settings)
        mpc.settings.reference_tracking = true
    else
        mpc.settings = settings
    end

    range = ParameterRange(mpc)
    range.xmax[:] .= 10 * ones(nx)
    range.xmin[:] .= -10 * ones(nx)
    range.rmax[:] .= 10 * ones(nx)
    range.rmin[:] .= -10 * ones(nx)

    x0 = zeros(nx)
    x0[1] = 3.0
    scenarios = [Scenario(x0; N = 15, r = zeros(nx))]
    return _finalize_example(:chained, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:mass_spring}, Np, Nc; settings=nothing, kwargs...)
    κ = get(kwargs, :κ, 1.0)
    λ = get(kwargs, :λ, 0.0)
    nm_kw = get(kwargs, :nm, nothing)
    nx_kw = get(kwargs, :nx, nothing)
    nm = if !isnothing(nm_kw)
        Int(nm_kw)
    elseif isnothing(nx_kw)
        1
    else
        nx = Int(nx_kw)
        max(1, fld(iseven(nx) ? nx : nx - 1, 2))
    end
    nx = 2 * nm

    Fx = diagm(1 => κ * ones(nm - 1), -1 => κ * ones(nm - 1), 0 => -2κ * ones(nm))
    Fv = diagm(1 => λ * ones(nm - 1), -1 => λ * ones(nm - 1), 0 => -2λ * ones(nm))
    A = [zeros(nm, nm) Matrix(I, nm, nm);
         Fx Fv]
    B = [zeros(nm, 1);
         1;
         zeros(nm - 1, 1)]
    C = Matrix(I, 2 * nm, 2 * nm)
    Ts = 0.5
    F, G = zoh(A, B, Ts)

    mpc = MPC(F, G; C, Np, Nc, Ts)
    set_objective!(mpc; Q = 100 * ones(nx), R = [1.0], Rr = [0.0])
    set_bounds!(mpc, umin = [-0.5], umax = [0.5])
    add_constraint!(mpc, Ax = Matrix(I, nm, 2 * nm), lb = -4 * ones(nm), ub = 4 * ones(nm), ks = 2:mpc.Nc)

    if isnothing(settings)
        mpc.settings.reference_tracking = false
    else
        mpc.settings = settings
    end

    range = ParameterRange(mpc)
    range.xmax[:] .= 4 * ones(nx)
    range.xmin[:] .= -4 * ones(nx)

    x0 = zeros(nx)
    x0[1] = 1.0
    scenarios = [Scenario(x0; N = 15)]
    return _finalize_example(:mass_spring, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:nonlinear_demo}, Np, Nc; settings=nothing, kwargs...)
    Ts = 0.2
    F = [0.8187 zeros(1, 4);
         0.1474 0.6550 -0.1637 0.0489 0.4878;
         0.01637 0.1637 0.9825 3.43e-3 0.0523;
         zeros(1, 3) 0.8013 -0.1801;
         zeros(1, 3) 0.1801 0.9813]
    G = [0.1813 0 0;
         0.0163 0.1637 3.43e-3;
         1.14e-3 0.0175 1.77e-4;
         0 0 0.1801;
         0 0 0.0186]
    C = [1.0 0 0 0 0; 0 1 2 0 0]

    mpc = MPC(F, G; C, Np, Nc, Ts)
    set_objective!(mpc; Q = [1.0, 1.0], R = zeros(3), Rr = (1e-1 * [1, 1, 1]).^2)
    set_bounds!(mpc, umin = [-3.0, 2, 2], umax = [3.0, 2, 2])

    if isnothing(settings)
        mpc.settings.reference_tracking = true
    else
        mpc.settings = settings
    end

    range = ParameterRange(mpc)
    range.xmax[:] .= [2; ones(4)]
    range.xmin[:] .= -0.5 * ones(5)
    range.rmax[:] .= 10 * ones(2)
    range.rmin[:] .= -10 * ones(2)

    scenarios = [Scenario([0.5, 0.0, 0.0, 0.0, 0.0]; N = 15, r = [1.0, 0.0])]
    return _finalize_example(:nonlinear_demo, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:invpend_contact}, Np, Nc; settings=nothing, kwargs...)
    nwalls = min(Int(get(kwargs, :nwalls, 2)), 2)
    mc = get(kwargs, :mc, 1.0)
    mp = get(kwargs, :mp, 1.0)
    l = get(kwargs, :l, 1.0)
    d = get(kwargs, :d, 0.5)
    g = 10.0
    κ = get(kwargs, :κ, 100.0)
    ν = get(kwargs, :ν, get(kwargs, :v, 10.0))

    A = [0 0 1 0;
         0 0 0 1;
         0 (mp*g/mc) 0 0;
         0 (mc + mp)*g/(mc*l) 0 0]
    B = [0 0 0;
         0 0 0;
         1 / mc 0 0;
         1 / (mc * l) -1 / (mp * l) 1 / (mp * l)]

    C = Matrix{Float64}(I, 4, 4)
    Ts = 0.05
    F, G = zoh(A, B, Ts)

    uby = [d; pi / 10; 1; 1]
    lby = -uby
    δ2l, δ2u = -uby[1] + l * lby[2] - d, -lby[1] + l * uby[2] - d
    dotδ2l, dotδ2u = -uby[3] + l * lby[4], -lby[3] + l * uby[4]
    δ3l, δ3u = lby[1] - l * uby[2] - d, uby[1] - l * lby[2] - d
    dotδ3l, dotδ3u = lby[3] - l * uby[4], uby[3] - l * lby[4]

    u2l, u2u = κ * δ2l + ν * dotδ2l, κ * δ2u + ν * dotδ2u
    u3l, u3u = κ * δ3l + ν * dotδ3l, κ * δ3u + ν * dotδ3u

    ndelta = 4
    mld = MLDModel(F, G, zeros(4, ndelta), zeros(4, 0); C)
    mpc = MPC(mld; Np, Nc)
    Q = [1.0, 1, 1, 1]
    R = [1.0; 1e-4 * ones(6)]
    Rr = zeros(length(R))
    Qf, ~ = ared(mpc.model.F, mpc.model.G[:, 1], mpc.weights.R[1:1, 1:1], mpc.weights.Q)
    set_objective!(mpc; Q, R, Rr, Qf)
    set_input_bounds!(mpc, umin = [-1.0; 0; zeros(4)], umax = [1.0; 1e30; 1e30; ones(4)])

    if isnothing(settings)
        mpc.settings.reference_tracking = false
    else
        mpc.settings = settings
    end

    set_output_bounds!(mpc, ymin = lby, ymax = uby, ks = 2:mpc.Nc)

    function add_wall_contact_constraints!(mpc, uidx, gap_delta, force_delta,
                                           gap_ax, gap_c, gap_l, gap_u, force_ax, force_c, force_l, force_u)
        ks = 2:mpc.Nc

        add_indicator_constraint!(mpc, gap_delta;
                                  Ax = reshape(gap_ax, 1, :),
                                  c = [gap_c],
                                  m = gap_l,
                                  M = gap_u,
                                  ϵ = 0.0,
                                  sense = :ge,
                                  ks)
        add_indicator_constraint!(mpc, force_delta;
                                  Ax = reshape(force_ax, 1, :),
                                  c = [force_c],
                                  m = force_l,
                                  M = force_u,
                                  ϵ = 0.0,
                                  sense = :ge,
                                  ks)

        Ax_res = [zeros(2, mpc.model.nx);
                  -reshape(force_ax, 1, :);
                   reshape(force_ax, 1, :)]
        Au_res = zeros(4, 3)
        Au_res[:, uidx] .= [1.0, 1.0, 1.0, -1.0]
        Adelta_res = zeros(4, ndelta)
        Adelta_res[1, gap_delta] = -force_u
        Adelta_res[2, force_delta] = -force_u
        Adelta_res[3, force_delta] = -force_l
        Adelta_res[4, gap_delta] = force_u
        ub_res = [0.0, 0.0, -force_l - force_c, force_u - force_c]
        add_mld_constraint!(mpc; Ax = Ax_res, Au = Au_res, Adelta = Adelta_res, ub = ub_res, ks)
    end

    add_wall_contact_constraints!(mpc, 2, 1, 3,
                                  [-1, l, 0, 0], -d, δ2l, δ2u,
                                  [-κ, κ*l, -ν, ν*l], -κ * d, u2l, u2u)
    if nwalls == 2
        add_wall_contact_constraints!(mpc, 3, 2, 4,
                                      [1, -l, 0, 0], -d, δ3l, δ3u,
                                      [κ, -κ*l, ν, -ν*l], -κ * d, u3l, u3u)
    end

    range = ParameterRange(mpc)
    range.xmax[:] .= 20 * ones(4)
    range.xmin[:] .= -20 * ones(4)

    scenarios = [Scenario([0.0, 0.05, 0.0, 0.0]; N = 20)]
    return _finalize_example(:invpend_contact, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:ballplate}, Np, Nc; settings=nothing, kwargs...)
    A = [0 1.0 0 0;
         0 0 700 0;
         0 0 0 1;
         0 0 0 -34.69]
    B = [0; 0; 0; 3.1119;;]
    Ts = 0.03
    C = [1.0 0 0 0]

    F, G = zoh(A, B, Ts)
    mpc = MPC(F, G; Ts, Np, Nc, C)
    set_bounds!(mpc, umin = [-10.0], umax = [10.0])
    xbounds = [30; 15; 15 * pi / 180; 1]
    add_constraint!(mpc; Ax = Matrix(I, 4, 4), lb = -xbounds, ub = xbounds, soft = false)
    set_objective!(mpc; Q = [100.0], R = [0.1], Rr = [0.0], Qf = [1.0])
    isnothing(settings) || (mpc.settings = settings)

    range = ParameterRange(mpc)
    range.xmax[:] = xbounds
    range.xmin[:] = -xbounds

    scenarios = [Scenario([10.0, 0.0, 0.0, 0.0]; T = 2.0, r = [0.0])]
    return _finalize_example(:ballplate, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:crazyflie}, Np, Nc; settings=nothing, kwargs...)
    mass = get(kwargs, :mass, 0.035)
    arm_length= get(kwargs, :arm_length, 0.046/1.414213562)
    J = get(kwargs, :J,[1.66e-5 0.83e-6 0.72e-6; 
                        0.83e-6 1.66e-5 1.8e-6; 
                        0.72e-6 1.8e-6 2.93e-5])
    thrustToTorque = get(kwargs,:thrustToTorque, 0.0008)
    kt = get(kwargs,:kt, 2.245365e-6*65536) #thrust coefficient (PWM scale: 2^16-1)

    diagonal_intertia = get(kwargs,:diagonal_intertia,true)
    diagonal_intertia && (J=Diagonal(J))

    g = 9.81
    km = kt*thrustToTorque # moment coefficient

    function cf_dynamics(x,u,d)
        # auxiliary functions
        hat(v) = [0 -v[3] v[2];v[3] 0 -v[1];-v[2] v[1] 0]
        L(q) = [q[1] -q[2:4]'; q[2:4] q[1]*I+hat(q[2:4])]
        T = cat(I(1),-I(3),dims=(1,2))
        H = [zeros(1,3); I]
        qtoQ(q) = H'*T*L(q)*T*L(q)*H

        # Extract states
        r = x[1:3]
        q = x[4:6]
        v = x[7:9]
        ω = x[10:12]

        qe = [1-q'*q;q] # true quaternion
        Q = qtoQ(qe)

        ṙ = v
        q̇ = 0.5*L(qe)*H*ω

        v̇ = [0; 0; -g] + (1/mass)*Q*[zeros(2,4); kt*ones(1,4)]*u 

        Cu = [(arm_length*kt)*[-1 -1 1 1; -1 1 1 -1]; km*[-1 1 -1 1]];
        ω̇ = J\(-hat(ω)*J*ω + Cu*u)

        return [ṙ; q̇[2:4]; v̇; ω̇]
    end

    # Steady state when hovering  
    x0 = zeros(12) 
    u0 = (mass*g/kt/4)*ones(4)

    # Start setting up MPC
    Ts = get(kwargs,:Ts,1/500)
    model = LinearMPC.Model(cf_dynamics,(x,u,d)->x,x0,u0,Ts)

    mpc = LinearMPC.MPC(model;Np,Nc);
    mpc.settings.reference_tracking = false;

    Q =[156.25, 156.25, 400,  # position
        2.777778, 2.777778, 1111.11111,  #angle
        4,4,4, # velocity
        4,4,25] # angular velocity
    R = 50*[1,1,1,1]
    set_objective!(mpc;Q=Q,R=R, Rr=0);

    set_bounds!(mpc,umin = zeros(4), umax=ones(4));
    set_terminal_cost!(mpc);
    set_prestabilizing_feedback!(mpc);

    if isnothing(settings)
        mpc.settings.reference_tracking = false
    else
        mpc.settings = settings
    end

    range = ParameterRange(mpc)
    range.xmax[:] .= 1
    range.xmin[:] .=-1

    x0 = zeros(12)
    x0[4] = 0.1
    x0[5] = -0.1
    scenarios = [Scenario(x0; T = 1.5)]
    return _finalize_example(:crazyflie, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:satellite}, Np, Nc; settings=nothing, kwargs...)
    A = [0.0 1 0; 0 0 0; 0 0 0]
    B = [0 0 0; 2.5 1 1; -10 0 0]

    mpc = MPC(A, B, 0.1; Np, Nc)
    set_objective!(mpc; Q = [0.5e4, 1e-2, 1e-1], R = [10, 10, 10], Rr = 0)
    set_bounds!(mpc; umin = [-Inf; 0; -1], umax = [Inf; 1; 0])
    set_binary_controls!(mpc, [2, 3])
    isnothing(settings) || (mpc.settings = settings)

    range = ParameterRange(mpc)
    scenarios = [Scenario(zeros(3); N = 20, r = [0.5, 0.0, 0.0])]
    return _finalize_example(:satellite, mpc, range; scenarios)
end

function _build_mpc_example(::Val{:rocket}, Np, Nc; settings=nothing, kwargs...)
    mass = get(kwargs,:mass, 530.406)
    inertia = get(kwargs,:intertia, 1209.5)
    l1 = get(kwargs,:l1, 2.8467)
    l2 = get(kwargs,:l2, 2.135)
    gravity = 9.81

    Ts = get(kwargs,:Ts, 0.01)

    scale, lander_scaling = 30.0, 4.0
    main_engine_thrust = (6.8e6 * lander_scaling^3) / scale^3
    side_engine_thrust = main_engine_thrust / 50.0
    max_nozzle_angle = 15.0 * pi / 180.0
    function rocket_dynamics(x,u,d)
        _, _, dx, dy, theta, dtheta = x
        u_main, u_side, u_gimbal = u

        main_force = clamp(u_main, 0.0, 1.0) * main_engine_thrust
        side_force = clamp(u_side, -1.0, 1.0) * side_engine_thrust
        phi = clamp(u_gimbal, -1.0, 1.0) * max_nozzle_angle

        ddx = (-main_force * sin(theta + phi) + side_force * cos(theta)) / mass
        ddy = (main_force * cos(theta + phi) + side_force * sin(theta)) / mass - gravity
        ddtheta = (-l1 * main_force * sin(phi) - l2 * side_force) / inertia
        return [dx, dy, ddx, ddy, dtheta, ddtheta]
    end

    x0 = zeros(6)
    u0 = [(mass * gravity) / main_engine_thrust, 0, 0]

    model = LinearMPC.Model(rocket_dynamics,(x,u,d)->x,x0,u0,Ts)
    mpc = LinearMPC.MPC(model;Np,Nc);
    set_objective!(mpc;Q=[2,1,5,5,1,5],R=0.1*ones(3), Rr=0);
    set_bounds!(mpc,umin = [0.0, -1.0,-1.0], umax=ones(3));

    if isnothing(settings)
        mpc.settings.reference_tracking = false
    else
        mpc.settings = settings
    end

    range = ParameterRange(mpc)
    scenarios = [Scenario([0.7, 1, 0.65, -15, 0.12, 0.05]; N = 250)]
    return _finalize_example(:rocket, mpc, range; scenarios)
end

function mpc_example(s::Union{AbstractString,Symbol}, Np, Nc = Np; params = Dict(), settings = nothing, kwargs...)
    id = _resolve_example_id(s)
    merged_kwargs = _merge_example_kwargs(params, kwargs)
    return _build_mpc_example(Val(id), Np, Nc; settings, merged_kwargs...)
end

function mpc_example(s::Union{AbstractString,Symbol}; params = Dict(), settings = nothing, kwargs...)
    id = _resolve_example_id(s)
    Np, Nc = _default_horizons(id)
    return mpc_example(id, Np, Nc; params, settings, kwargs...)
end

function mpc_examples(args...; kwargs...)
    example = mpc_example(args...; kwargs...)
    return example.mpc, example.range
end

function Simulation(example::MPCExample, scenario::Scenario)
    return Simulation(example.mpc, scenario)
end

function Simulation(example::MPCExample, scenario_id::Integer)
    1 <= scenario_id <= length(example.scenarios) || throw(BoundsError(example.scenarios, scenario_id))
    return Simulation(example, example.scenarios[scenario_id])
end

function Simulation(example::MPCExample; kwargs...)
    if !isempty(example.scenarios)
        return Simulation(example,1)
    else
        return Simulation(example.mpc; kwargs...)
    end
end
