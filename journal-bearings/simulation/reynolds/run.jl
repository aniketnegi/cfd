# Named-case CLI for the verified *cylindrical, unported, signed full-film* model.
include(joinpath(@__DIR__,"src/convergence.jl"))
include(joinpath(@__DIR__, "..", "..", "experiment.jl"))

using ArgParse, Printf

function run_settings()
    settings = ArgParseSettings(
        prog="julia --project=. simulation/reynolds/run.jl",
        description="Solve signed full-film Reynolds pressure on a named unported cylinder.",
        epilog="""New runs require --rpm and --viscosity-pa-s. Existing runs reuse
        saved inputs; --replace permits rerunning but not changing those inputs.
        Results are written under out/cases/EXPERIMENT/runs/RUN/.
        Negative pressure is retained; this model does not include cavitation.

        Example: --experiment demo-cylinder --run vg32 --rpm 2000 --viscosity-pa-s 0.02752

        Exit status: 0 success, 1 execution failure, 2 invalid input.""",
        exit_after_help=false,
        exc_handler=(_, err) -> throw(ArgumentError(err.text)))
    @add_arg_table! settings begin
        "--experiment"
            help = "existing geometry case"
            metavar = "NAME"
            arg_type = String
            required = true
            action = :append_arg
        "--run"
            help = "operating-condition run name"
            metavar = "NAME"
            arg_type = String
            required = true
            action = :append_arg
        "--rpm"
            help = "signed journal speed in revolutions per minute"
            metavar = "RPM"
            arg_type = Float64
            range_tester = isfinite
            action = :append_arg
        "--viscosity-pa-s"
            help = "positive dynamic viscosity in Pa s"
            metavar = "PA_S"
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
            action = :append_arg
        "--sweep"
            help = "also write a grid-convergence table"
            action = :count_invocations
        "--replace"
            help = "rerun with the same saved inputs"
            action = :count_invocations
    end
    add_arg_group!(settings, "output", "output"; exclusive=true)
    @add_arg_table! settings begin
        "--verbose", "-v"
            help = "show timings and full numerical diagnostics"
            action = :store_true
        "--quiet", "-q"
            help = "suppress routine output"
            action = :store_true
    end
    return settings
end

function parse_run_args(args; cases_root=DEFAULT_CASES_ROOT)
    parsed = parse_args(args, run_settings())
    isnothing(parsed) && return nothing # --help
    # Preserve rejection of repeated scientific inputs rather than silently
    # choosing the last occurrence.
    opts = Dict{String,Any}()
    for key in ("experiment", "run", "rpm", "viscosity-pa-s")
        values = parsed[key]
        length(values) <= 1 || throw(ArgumentError("repeated option: --$key"))
        isempty(values) || (opts["--$key"] = only(values))
    end
    for key in ("replace", "sweep")
        parsed[key] <= 1 || throw(ArgumentError("repeated --$key"))
    end
    replace, sweep = parsed["replace"] == 1, parsed["sweep"] == 1
    verbose, quiet = parsed["verbose"], parsed["quiet"]
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
    rpm, mu = get(opts, "--rpm", nothing), get(opts, "--viscosity-pa-s", nothing)
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
            existing, sweep, replace, verbose, quiet, experiment_sha256)
end

function main_run(args=ARGS; cases_root=DEFAULT_CASES_ROOT)
    cfg = parse_run_args(args; cases_root)
    isnothing(cfg) && return nothing
    execute_run(cfg)
end

function execute_run(cfg)
    if !cfg.quiet
        @printf(stderr, "run %s/%s grid=%dx%d rpm=%g viscosity=%g Pa s\n",
            cfg.exp.data["name"], cfg.run_name, cfg.exp.ntheta, cfg.exp.nz, cfg.rpm, cfg.mu)
        flush(stderr)
    end
    cfg.verbose && println(stderr, "mesh ", abspath(cfg.exp.mesh))
    solve_seconds = @elapsed r = solve_cylindrical_full_film(cfg.exp.p, cfg.exp.ntheta, cfg.exp.nz;
                                    rpm=cfg.rpm, viscosity_pa_s=cfg.mu)
    cfg.verbose && @printf(stderr, "solve %.3f s\n", solve_seconds)
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
        if !cfg.quiet
            println(stderr, "sweep ", join(["$(n)x$(z)" for (n,z) in reverse(levels)], ", "))
            flush(stderr)
        end
        tmp_sweep, io_sweep = mktemp(cfg.run_dir)
        close(io_sweep)
        try
            sweep_seconds = @elapsed write_convergence(tmp_sweep, cfg.exp.p, reverse(levels);
                              rpm=cfg.rpm, viscosity_pa_s=cfg.mu)
            chmod(tmp_sweep, stat(cfg.config_path).mode & 0o777)
            mv(tmp_sweep, cfg.sweep_path; force=cfg.replace && ispath(cfg.sweep_path))
            cfg.verbose && @printf(stderr, "sweep %.3f s\n", sweep_seconds)
        finally
            ispath(tmp_sweep) && rm(tmp_sweep)
        end
    end
    if !cfg.quiet
        @printf(stderr, "pressure min=%.6g max=%.6g Pa (signed)\n", extrema(r.pressure_pa)...)
        @printf(stderr, "residual relative_l1=%.3g\n", r.diagnostics.relative_l1_cell_imbalance)
        println(stderr, "output ", abspath(cfg.run_dir))
    end
    if cfg.verbose
        for key in sort!(collect(keys(r.diagnostics)); by=string)
            println(stderr, key, "=", getproperty(r.diagnostics, key))
        end
    end
    return nothing
end

function run_cli(args=ARGS; cases_root=DEFAULT_CASES_ROOT)
    cfg = try
        parse_run_args(args; cases_root)
    catch err
        err isa InterruptException && return 130
        println(stderr, "run.jl: ", sprint(showerror, err))
        if err isa ArgumentError || err isa TOML.ParserError
            println(stderr, "Try --help for usage.")
            return 2
        end
        return 1
    end
    isnothing(cfg) && return 0
    try
        execute_run(cfg)
    catch err
        if err isa InterruptException
            println(stderr, "run.jl: interrupted")
            return 130
        end
        println(stderr, "run.jl: ", sprint(showerror, err))
        if cfg.verbose
            Base.show_backtrace(stderr, catch_backtrace())
            println(stderr)
        end
        return 1
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(run_cli())
end
