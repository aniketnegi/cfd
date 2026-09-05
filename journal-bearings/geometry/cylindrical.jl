using Gmsh: Gmsh, gmsh


struct BearingParameters
    journal_radius::Float64
    radial_clearance::Float64
    length::Float64
    eccentricity_ratio::Float64

    housing_radius::Float64
    eccentricity::Float64
    journal_offset::NTuple{3,Float64}

    function BearingParameters(;
        journal_radius,
        radial_clearance,
        length,
        eccentricity_ratio,
    )

        jr = Float64(journal_radius)
        rc = Float64(radial_clearance)
        l = Float64(length)
        er = Float64(eccentricity_ratio)

        isfinite(jr) && jr > 0 ||
            throw(ArgumentError("journal_radius must be positive and finite, got $jr"))
        isfinite(rc) && rc > 0 ||
            throw(ArgumentError("radial_clearance must be positive and finite, got $rc"))
        isfinite(l) && l > 0 ||
            throw(ArgumentError("length must be positive and finite, got $l"))
        isfinite(er) && 0 <= er < 1 ||
            throw(ArgumentError("eccentricity_ratio must be in [0, 1), got $er"))

        hr = jr + rc
        ecc = rc * er
        offset = (ecc, 0.0, 0.0)

        return new(jr, rc, l, er, hr, ecc, offset)
    end
end

function makebearing(params::BearingParameters)
    housing_tag =
        gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, params.length, params.housing_radius)

    journal_tag = gmsh.model.occ.addCylinder(
        params.journal_offset...,
        0,
        0,
        params.length,
        params.journal_radius,
    )

    out_dimtags, _ = gmsh.model.occ.cut([(3, housing_tag)], [(3, journal_tag)])

    length(out_dimtags) == 1 ||
        error("Expected one Boolean result, got $(length(out_dimtags))")

    dimension, fluid_tag = only(out_dimtags)
    dimension == 3 || error("Expected a volume, got dimension $dimension")

    gmsh.model.occ.synchronize()

    # Volume Check
    volume = gmsh.model.occ.getMass(3, fluid_tag)
    expected_volume =
        π * (params.housing_radius^2 - params.journal_radius^2) * params.length

    isapprox(volume, expected_volume; rtol = 1e-9, atol = 0.0) ||
        error("Volume mismatch: expected $expected_volume mm³, got $volume mm³")


    # Bounds Check
    bounds = gmsh.model.getBoundingBox(3, fluid_tag)

    expected_bounds = (
        -params.housing_radius,
        -params.housing_radius,
        0.0,
        params.housing_radius,
        params.housing_radius,
        params.length,
    )

    all(isapprox.(bounds, expected_bounds; atol = 2e-7, rtol = 0.0)) ||
        error("Bounds mismatch: expected $expected_bounds, got $bounds")

    return fluid_tag
end


function main()
    params = BearingParameters(
        journal_radius = 50.0,
        radial_clearance = 0.05,
        length = 60.0,
        eccentricity_ratio = 0.6,
    )

    output_path = normpath(joinpath(@__DIR__, "..", "out", "bearing.brep"))

    Gmsh.initialize()

    try
        gmsh.model.add("bearing")

        makebearing(params)

        mkpath(dirname(output_path))
        gmsh.write(output_path)

        isfile(output_path) || error("Gmsh did not create $output_path")

        filesize(output_path) > 0 || error("Gmsh created an empty BREP at $output_path")

        println("Saved BREP to $output_path")
    finally
        Gmsh.finalize()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
