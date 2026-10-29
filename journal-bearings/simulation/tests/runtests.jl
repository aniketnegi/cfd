# Run each suite in a fresh process: includes define shared types and globals.
const SUITES = [
    "full_film_checks.jl", "fitted_full_film_checks.jl", "swift_stieber_checks.jl",
    "fitted_swift_stieber_checks.jl", "pressure_sampling_checks.jl", "checks.jl",
    "conical_checks.jl", "jfo_slider_checks.jl", "conical_jfo_checks.jl",
    "cylindrical_jfo_checks.jl", "xi2024_jfo_checks.jl", "fitted_jfo_checks.jl",
    "elasticity_checks.jl", "fsi_checks.jl",
    "journal_checks.jl", "both_fsi_checks.jl",
    "conical_equilibrium_checks.jl", "run_checks.jl",
]

function main_tests(args=ARGS)
    selected = String[]
    verbose = quiet = false
    options = true
    for arg in args
        if options && arg == "--"
            options = false
        elseif options && arg in ("-h", "--help")
            println("Usage: julia --project=. simulation/tests/runtests.jl [OPTIONS] [SUITE.jl ...]")
            println("Run all suites by default, each in a separate process.")
            println("  -v, --verbose  show raw Julia test output")
            println("  -q, --quiet    show only failures")
            println("  -h, --help     show this help")
            println("Suites:")
            foreach(suite -> println("  ", suite), SUITES)
            return 0
        elseif options && arg in ("-v", "--verbose")
            verbose = true
        elseif options && arg in ("-q", "--quiet")
            quiet = true
        elseif arg in SUITES
            push!(selected, arg)
        else
            println(stderr, "runtests: unknown argument: ", arg)
            return 2
        end
    end
    if verbose && quiet
        println(stderr, "runtests: --verbose and --quiet cannot be combined")
        return 2
    end
    isempty(selected) && append!(selected, SUITES)
    root = normpath(joinpath(@__DIR__, "../.."))
    failed = String[]
    for suite in selected
        progress = stderr isa Base.TTY && !quiet && !verbose
        if progress
            print(stderr, "run ", suite, " ")
        elseif !quiet
            println(stderr, "run ", suite)
        end
        flush(stderr)
        command = `$(Base.julia_cmd()) --color=no --project=$root --check-bounds=yes $(joinpath(@__DIR__, suite))`
        passed = mktemp() do path, log
            started = time_ns()
            heartbeat = progress ? Timer(1; interval=1) do timer
                print(stderr, '.'); flush(stderr)
            end : nothing
            process = try
                if verbose
                    run(pipeline(ignorestatus(command); stdout=stderr, stderr=stderr))
                else
                    run(pipeline(ignorestatus(command); stdout=log, stderr=log))
                end
            finally
                isnothing(heartbeat) || close(heartbeat)
            end
            elapsed = round(Int, (time_ns() - started) / 1e9)
            ok = success(process)
            if progress
                println(stderr, ok ? " [OK] [$(elapsed)s]" : " [FAIL] [$(elapsed)s]")
                flush(stderr)
            end
            if !ok
                println(stderr, "FAIL ", suite, " (exit ", process.exitcode, ")")
                if !verbose
                    flush(log)
                    seekstart(log)
                    for line in eachline(log)
                        println(stderr, line)
                    end
                end
            elseif !quiet && !progress
                println(stderr, "ok ", suite)
            end
            ok
        end
        passed || push!(failed, suite)
    end
    if !isempty(failed)
        println(stderr, "runtests: ", length(failed), " suite(s) failed")
        return 1
    end
    quiet || println(stderr, "runtests: ", length(selected), " suite(s) passed")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_tests())
end
