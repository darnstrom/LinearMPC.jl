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
    B = [B zeros(4, 4)]

    C, D = Matrix{Float64}(I, 4, 4), zeros(4, 1)
    Ts = 0.05
    F, G = zoh(A, B, Ts)

    mpc = MPC(F, G; C, Np, Nc)
    Q = [1.0, 1, 1, 1]
    R = [1.0; 1e-4 * ones(6)]
    Rr = zeros(7)
    Qf, ~ = ared(mpc.model.F, mpc.model.G[:, 1], mpc.weights.R[1:1, 1:1], mpc.weights.Q)
    set_objective!(mpc; Q, R, Rr, Qf)
    set_bounds!(mpc, umin = [-1.0; 0; zeros(4)], umax = [1.0; 1e30; 1e30; ones(4)])
    set_binary_controls!(mpc, collect(4:7))

    if isnothing(settings)
        mpc.settings.reference_tracking = false
    else
        mpc.settings = settings
    end

    uby = [d; pi / 10; 1; 1]
    lby = -uby
    set_output_bounds!(mpc, ymin = lby, ymax = uby, ks = 2:mpc.Nc)

    δ2l, δ2u = -uby[1] + l * lby[2] - d, -lby[1] + l * uby[2] - d
    dotδ2l, dotδ2u = -uby[3] + l * lby[4], -lby[3] + l * uby[4]
    δ3l, δ3u = lby[1] - l * uby[2] - d, uby[1] - l * lby[2] - d
    dotδ3l, dotδ3u = lby[3] - l * uby[4], uby[3] - l * lby[4]

    u2l, u2u = κ * δ2l + ν * dotδ2l, κ * δ2u + ν * dotδ2u
    u3l, u3u = κ * δ3l + ν * dotδ3l, κ * δ3u + ν * dotδ3u

    Ax = [-1 l 0 0;
          1 -l 0 0;
          -κ κ*l -ν ν*l;
          κ -κ*l ν -ν*l;
          zeros(2, 4);
          κ -κ*l ν -ν*l;
          -κ κ*l -ν ν*l]
    Au2 = [0 0 0 -δ2u 0 0 0;
           0 0 0 -δ2l 0 0 0;
           0 0 0 0 0 -u2u 0;
           0 0 0 0 0 -u2l 0;
           0 1 0 -u2u 0 0 0;
           0 1 0 0 0 -u2u 0;
           0 1 0 0 0 -u2l 0;
           0 -1 0 u2u 0 0 0]
    Au3 = [0 0 0 0 -δ3u 0 0;
           0 0 0 0 -δ3l 0 0;
           0 0 0 0 0 0 -u3u;
           0 0 0 0 0 0 -u3l;
           0 0 1 0 -u3u 0 0;
           0 0 1 0 0 0 -u3u;
           0 0 1 0 0 0 -u3l;
           0 0 -1 0 u3u 0 0]
    bg2 = [d;
           -δ2l - d;
           κ * d;
           -κ * d - u2l;
           0;
           0;
           -u2l - κ * d;
           u2u + κ * d]
    bg3 = [d;
           -δ3l - d;
           κ * d;
           -κ * d - u3l;
           0;
           0;
           -u3l - κ * d;
           u3u + κ * d]

    add_constraint!(mpc, Au = Au2, Ax = Ax, ub = bg2, ks = 2:mpc.Nc)
    if nwalls == 2
        add_constraint!(mpc, Au = Au3, Ax = -Ax, ub = bg3, ks = 2:mpc.Nc)
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
    F = [1.0 0.0 0.0 0.0000000 0.0009810 0.0000000 0.0100000 0.0000000 0.0000000 0.0000000 0.0000016 0.0000000;
         0.0 1.0 0.0 -0.0009810 0.0000000 0.0000000 0.0000000 0.0100000 0.0000000 -0.0000016 0.0000000 0.0000000;
         0.0 0.0 1.0 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0100000 0.0000000 0.0000000 0.0000000;
         0.0 0.0 0.0 1.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0050000 0.0000000 0.0000000;
         0.0 0.0 0.0 0.0000000 1.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0050000 0.0000000;
         0.0 0.0 0.0 0.0000000 0.0000000 1.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0050000;
         0.0 0.0 0.0 0.0000000 0.1962000 0.0000000 1.0000000 0.0000000 0.0000000 0.0000000 0.0004905 0.0000000;
         0.0 0.0 0.0 -0.1962000 0.0000000 0.0000000 0.0000000 1.0000000 0.0000000 -0.0004905 0.0000000 0.0000000;
         0.0 0.0 0.0 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 1.0000000 0.0000000 0.0000000 0.0000000;
         0.0 0.0 0.0 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 1.0000000 0.0000000 0.0000000;
         0.0 0.0 0.0 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 1.0000000 0.0000000;
         0.0 0.0 0.0 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 1.0000000]
    G = [-0.0000011 0.0000012 0.0000011 -0.0000012;
         0.0000011 0.0000012 -0.0000011 -0.0000012;
         0.0002102 0.0002102 0.0002102 0.0002102;
         -0.0068839 -0.0075809 0.0068916 0.0075732;
         -0.0069177 0.0076070 0.0069392 -0.0076285;
         0.0004937 -0.0001806 -0.0006961 0.0003830;
         -0.0004524 0.0004975 0.0004538 -0.0004989;
         0.0004502 0.0004958 -0.0004507 -0.0004953;
         0.0420429 0.0420429 0.0420429 0.0420429;
         -2.7535461 -3.0323404 2.7566264 3.0292601;
         -2.7670702 3.0427842 2.7756950 -3.0514090;
         0.1974771 -0.0722364 -0.2784376 0.1531969]

    Q = 1 ./ ([0.1; 0.1; 0.1; 0.5; 0.5; 0.03; 0.5; 0.5; 0.5; 0.7; 0.7; 0.2].^2)
    R = 1 ./ (([0.5; 0.5; 0.5; 0.5] / 6).^2)
    Ts = 1e-2
    u0 = [0.5833333520642209, 0.5833333520642209, 0.5833333520642209, 0.5833333520642209]

    mpc = MPC(F, G; Ts, Np, Nc)
    set_bounds!(mpc, umin = -u0, umax = ones(4) - u0)
    set_objective!(mpc; Q, R, Rr = 0)

    if isnothing(settings)
        mpc.settings.reference_tracking = false
    else
        mpc.settings = settings
    end

    range = ParameterRange(mpc)
    range.xmax[:] .= 5
    range.xmin[:] .= -5

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

Simulation(example::MPCExample; kwargs...) = Simulation(example.mpc; kwargs...)
