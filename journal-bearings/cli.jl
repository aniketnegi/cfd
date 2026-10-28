"""Small command-line helpers; no numerical dependencies."""
module BearingCLI

function study_options(args; full=false, default_output=nothing, output=true,
                       description="Run a numerical verification study.")
    delimiter = findfirst(==("--"), args)
    option_args = isnothing(delimiter) ? args : args[1:delimiter-1]
    if any(arg -> arg in ("-h", "--help"), option_args)
        destination = output ? (isnothing(default_output) ? " DIRECTORY" : " [DIRECTORY]") : ""
        println("Usage: ", basename(PROGRAM_FILE), " [OPTIONS]", destination)
        println(description)
        output && println("DIRECTORY must be fresh; use -- before a path beginning with '-'.")
        isnothing(default_output) || println("Default directory: ", default_output)
        full && println("  --full         run the full study instead of representative cases")
        println("  -v, --verbose  show case diagnostics and detailed errors")
        println("  -q, --quiet    show only warnings and errors")
        println("  -h, --help     show help")
        println("Exit status: 0 success, 1 execution/check failure, 2 invalid arguments, 130 interrupt.")
        return nothing
    end
    directory = nothing
    verbose = quiet = all_cases = false
    options = true
    for arg in args
        if options && arg == "--"
            options = false
        elseif options && arg in ("-v", "--verbose")
            verbose = true
        elseif options && arg in ("-q", "--quiet")
            quiet = true
        elseif options && arg == "--full" && full
            all_cases && throw(ArgumentError("repeated --full"))
            all_cases = true
        elseif options && startswith(arg, "-")
            throw(ArgumentError("unknown option: $arg"))
        elseif output && isnothing(directory)
            directory = arg
        else
            throw(ArgumentError("unexpected argument: $arg"))
        end
    end
    verbose && quiet && throw(ArgumentError("--verbose and --quiet cannot be combined"))
    directory = isnothing(directory) ? default_output : directory
    if output
        (isnothing(directory) || isempty(directory)) && throw(ArgumentError("a fresh output directory is required"))
        (ispath(directory) || islink(directory)) && throw(ArgumentError("output already exists: $directory"))
    end
    return (; output=directory, full=all_cases, verbose, quiet)
end

function report(f; verbose=false, quiet=false, label=basename(PROGRAM_FILE))
    progress = stderr isa Base.TTY && !quiet && !verbose
    quiet || (progress ? print(stderr, "run $label ") : println(stderr, "run $label"))
    flush(stderr)
    started = time_ns()
    return mktemp() do path, log
        heartbeat = progress ? Timer(1; interval=1) do timer
            print(stderr, '.'); flush(stderr)
        end : nothing
        stop_progress() = isnothing(heartbeat) || close(heartbeat)
        try
            redirect_stdout(verbose ? stderr : log) do
                f()
            end
            stop_progress()
            elapsed = round((time_ns() - started) / 1e9; digits=1)
            quiet || println(stderr, progress ? " [OK] [$(elapsed)s]" : "ok $label [$(elapsed)s]")
            return 0
        catch err
            stop_progress()
            progress && println(stderr, " [FAIL]")
            if !verbose
                flush(log); seekstart(log)
                for line in eachline(log)
                    println(stderr, line)
                end
            end
            println(stderr, label, ": ", sprint(showerror, err))
            verbose && (Base.show_backtrace(stderr, catch_backtrace()); println(stderr))
            return err isa InterruptException ? 130 : 1
        finally
            stop_progress()
        end
    end
end

function study_main(f, args=ARGS; kwargs...)
    cfg = try
        study_options(args; kwargs...)
    catch err
        err isa InterruptException && return 130
        err isa ArgumentError || rethrow()
        println(stderr, basename(PROGRAM_FILE), ": ", sprint(showerror, err), ". Try --help.")
        return 2
    end
    isnothing(cfg) && return 0
    result = nothing
    code = report(; cfg.verbose, cfg.quiet) do
        result = f(cfg)
    end
    directory = isnothing(cfg.output) && result isa AbstractString ? result : cfg.output
    if code == 0 && !cfg.quiet && !isnothing(directory)
        println(stderr, "output ", abspath(directory))
    end
    return code
end

end
