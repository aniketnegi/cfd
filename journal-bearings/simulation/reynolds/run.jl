# Named-case CLI for the verified *cylindrical, unported, signed full-film* model.
include(joinpath(@__DIR__, "convergence.jl"))
include(joinpath(@__DIR__, "..", "..", "experiment.jl"))

const RUN_USAGE = """NAME
    run.jl - run cylindrical signed full-film Reynolds on a named experiment

SYNOPSIS
    julia --project=. simulation/reynolds/run.jl --experiment NAME --run NAME
        [--rpm RPM --viscosity-pa-s PA_S] [--sweep] [--replace]

DESCRIPTION
    Read geometry and resolution from out/cases/NAME/experiment.toml and
    verify the saved mesh checksum. Only an unported zero-angle cylinder is
    supported. New runs require rpm and dynamic viscosity (Pa s); later
    invocations may load them from runs/NAME/run.toml. An existing result is
    never overwritten without --replace, and a saved run's inputs cannot
    change under the same name. Positive/negative signed pressure is NOT JFO.
    Options accept --key VALUE or --key=VALUE.

REQUIRED OPTIONS
    --experiment=NAME          existing named geometry case
    --run=NAME                 safe identifier for operating conditions

OPTIONS
    --rpm=RPM                  signed journal speed, required for new runs
    --viscosity-pa-s=PA_S      positive dynamic viscosity, required for new runs
    --sweep                   also write a grid-convergence TSV for this run
    --replace                 explicitly rerun the identical saved inputs
    -h, --help                display this help and exit

EXIT STATUS
    0    Success.
    2    Invalid arguments.
"""

function parse_run_args(args; cases_root=DEFAULT_CASES_ROOT)
    opts = Dict{String,String}()
    replace = false
    sweep = false
    i = 1
    while i <= length(args)
        argument = args[i]
        if argument in ("--replace", "--sweep")
            if argument == "--replace"
                replace && throw(ArgumentError("repeated --replace"))
                replace = true
            else
                sweep && throw(ArgumentError("repeated --sweep"))
                sweep = true
            end
            i += 1
            continue
        end
        key, value = if startswith(argument, "--") && occursin('=', argument)
            i += 1
            split(argument, '='; limit=2)
        else
            i < length(args) && !startswith(args[i+1], "--") ||
                throw(ArgumentError("missing value for $argument"))
            i += 2
            (argument, args[i-1])
        end
        key in ("--experiment", "--run", "--rpm", "--viscosity-pa-s") ||
            throw(ArgumentError("unknown option: $key"))
        haskey(opts, key) && throw(ArgumentError("repeated option: $key"))
        isempty(value) && throw(ArgumentError("missing value for $key"))
        opts[key] = value
    end
    all(haskey(opts, key) for key in ("--experiment", "--run")) ||
        throw(ArgumentError("--experiment and --run are required"))
    exp = load_experiment(opts["--experiment"]; cases_root)
    exp.variant == "unported" || throw(ArgumentError("this runner supports unported cylinders only; surface-feed geometry does not define supply conditions"))
    iszero(exp.p.semicone_angle_deg) ||
        throw(ArgumentError("this named runner supports only zero-angle cylinders"))
    exp.ntheta >= 4 && exp.nz >= 2 || throw(ArgumentError("cylindrical solver requires grid ≥4×2"))
    sweep && !(iseven(exp.ntheta) && iseven(exp.nz) &&
               exp.ntheta ÷ 2 >= 4 && exp.nz ÷ 2 >= 2) &&
        throw(ArgumentError("--sweep requires at least two dyadic grids ≥4×2"))
    run_name = safe_case_name(opts["--run"])
    islink(joinpath(exp.dir, "runs")) &&
        throw(ArgumentError("runs directory must not be a symlink"))
    run_dir = joinpath(exp.dir, "runs", run_name)
    islink(run_dir) && throw(ArgumentError("run directory must not be a symlink: $run_dir"))
    config_path, result_path = joinpath(run_dir, "run.toml"), joinpath(run_dir, "full-film.vtk")
    sweep_path = joinpath(run_dir, "convergence.tsv")
    numeric(key) = begin
        haskey(opts, key) || return nothing
        v = tryparse(Float64, opts[key])
        !isnothing(v) && isfinite(v) || throw(ArgumentError("invalid finite number for $key"))
        return v
    end
    rpm, mu = numeric("--rpm"), numeric("--viscosity-pa-s")
    isnothing(mu) || mu > 0 || throw(ArgumentError("--viscosity-pa-s must be positive"))
    experiment_sha256 = mesh_hash(joinpath(exp.dir, "experiment.toml"))
    existing = isfile(config_path)
    if existing
        saved = TOML.parsefile(config_path)
        get(saved, "schema_version", nothing) === 1 &&
        get(saved, "experiment", nothing) == opts["--experiment"] &&
        get(saved, "name", nothing) == run_name &&
        get(saved, "model", nothing) == "cylindrical-signed-full-film" &&
        get(saved, "experiment_sha256", nothing) == experiment_sha256 &&
        get(saved, "axial_pressure_gauge_pa", nothing) == 0.0 &&
        get(saved, "circumferential_boundary", nothing) == "periodic" &&
        get(saved, "cavitation", nothing) == "none" ||
            throw(ArgumentError("invalid or stale run manifest: $config_path"))
        saved_rpm, saved_mu = get(saved, "rpm", nothing), get(saved, "viscosity_pa_s", nothing)
        saved_rpm isa Real && !(saved_rpm isa Bool) && isfinite(saved_rpm) &&
            saved_mu isa Real && !(saved_mu isa Bool) &&
            isfinite(saved_mu) && saved_mu > 0 || throw(ArgumentError("invalid run inputs: $config_path"))
        !isnothing(rpm) && rpm != saved_rpm &&
            throw(ArgumentError("rpm differs from saved run; choose a new --run name"))
        !isnothing(mu) && mu != saved_mu &&
            throw(ArgumentError("viscosity differs from saved run; choose a new --run name"))
        rpm, mu = Float64(saved_rpm), Float64(saved_mu)
    else
        isnothing(rpm) && throw(ArgumentError("new run requires --rpm"))
        isnothing(mu) && throw(ArgumentError("new run requires --viscosity-pa-s"))
        ispath(run_dir) && throw(ArgumentError("run directory exists without run.toml: $run_dir"))
        replace && throw(ArgumentError("--replace requires an existing run"))
    end
    ispath(result_path) && !replace &&
        throw(ArgumentError("result already exists: $result_path (use --replace to rerun)"))
    sweep && ispath(sweep_path) && !replace &&
        throw(ArgumentError("sweep already exists: $sweep_path (use --replace to rerun)"))
    return (; exp, run_name, run_dir, config_path, result_path, sweep_path, rpm, mu,
            existing, sweep, replace, experiment_sha256)
end

function main_run(args=ARGS; cases_root=DEFAULT_CASES_ROOT)
    if "--help" in args || "-h" in args
        print(RUN_USAGE)
        return
    end
    cfg = parse_run_args(args; cases_root)
    r = solve_cylindrical_full_film(cfg.exp.p, cfg.exp.ntheta, cfg.exp.nz;
                                    rpm=cfg.rpm, viscosity_pa_s=cfg.mu)
    if !cfg.existing
        mkpath(joinpath(cfg.exp.dir, "runs"))
        mkdir(cfg.run_dir)
        data = Dict("schema_version" => 1, "experiment" => cfg.exp.data["name"],
                    "name" => cfg.run_name, "model" => "cylindrical-signed-full-film",
                    "experiment_sha256" => cfg.experiment_sha256,
                    "rpm" => cfg.rpm, "viscosity_pa_s" => cfg.mu,
                    "axial_pressure_gauge_pa" => 0.0,
                    "circumferential_boundary" => "periodic", "cavitation" => "none")
        open(cfg.config_path, "w") do io
            TOML.print(io, data; sorted=true)
        end
    end
    # Write to a sibling temporary file, then rename: --replace must not
    # destroy a complete previous VTK if the new export fails midway.
    tmp, io = mktemp(cfg.run_dir)
    close(io)
    try
        write_pressure_vtk(tmp, r)
        chmod(tmp, stat(cfg.exp.mesh).mode & 0o777)
        mv(tmp, cfg.result_path; force=cfg.existing && ispath(cfg.result_path))
    finally
        ispath(tmp) && rm(tmp)
    end
    diagnostics_path = joinpath(cfg.run_dir, "diagnostics.toml")
    data = Dict{String,Any}(string(k) => (v isa Tuple ? collect(v) : v)
        for (k,v) in pairs(r.diagnostics))
    data["min_pressure_gauge_pa"] = minimum(r.pressure_pa)
    data["max_pressure_gauge_pa"] = maximum(r.pressure_pa)
    open(diagnostics_path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    if cfg.sweep
        levels = [(cfg.exp.ntheta, cfg.exp.nz)]
        nt, nz = first(levels)
        while length(levels) < 5 && iseven(nt) && iseven(nz) && nt ÷ 2 >= 4 && nz ÷ 2 >= 2
            nt, nz = nt ÷ 2, nz ÷ 2
            push!(levels, (nt, nz))
        end
        length(levels) >= 2 || throw(ArgumentError("--sweep requires at least two dyadic grids ≥4×2"))
        tmp_sweep, io_sweep = mktemp(cfg.run_dir)
        close(io_sweep)
        try
            write_convergence(tmp_sweep, cfg.exp.p, reverse(levels);
                              rpm=cfg.rpm, viscosity_pa_s=cfg.mu)
            chmod(tmp_sweep, stat(cfg.config_path).mode & 0o777)
            mv(tmp_sweep, cfg.sweep_path; force=cfg.replace && ispath(cfg.sweep_path))
            println("Wrote $(abspath(cfg.sweep_path))")
        finally
            ispath(tmp_sweep) && rm(tmp_sweep)
        end
    end
    println("Wrote $(abspath(cfg.result_path)) and $(abspath(diagnostics_path))")
    println("Signed full-film pressure [Pa]: $(extrema(r.pressure_pa)); negative values are not cavitation")
    println("Conservation diagnostics: $(r.diagnostics)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main_run()
    catch err
        if err isa ArgumentError
            println(stderr, "run.jl: $(sprint(showerror, err))")
            println(stderr, "Try 'julia --project=. simulation/reynolds/run.jl --help' for more information.")
            exit(2)
        end
        rethrow()
    end
end
