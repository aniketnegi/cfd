struct BearingParameters
    journal_midpoint_radius::Float64
    semicone_angle_deg::Float64
    radial_clearance::Float64
    bearing_length::Float64
    eccentricity_ratio::Float64

    housing_midpoint_radius::Float64
    eccentricity::Float64
    journal_offset::NTuple{3,Float64}
    minimum_gap::Float64
    maximum_gap::Float64


    # for perf reasons, calculating and caching here
    slope::Float64
    half_bearing_length::Float64

    function BearingParameters(;
        journal_midpoint_radius,
        semicone_angle_deg = 0.0,
        radial_clearance,
        bearing_length,
        eccentricity_ratio,
    )
        jmr = Float64(journal_midpoint_radius)
        rc = Float64(radial_clearance)
        bl = Float64(bearing_length)
        er = Float64(eccentricity_ratio)
        α = Float64(semicone_angle_deg)

        isfinite(jmr) && jmr > 0 || throw(
            ArgumentError("journal_midpoint_radius must be positive and finite, got $jmr"),
        )
        isfinite(rc) && rc > 0 ||
            throw(ArgumentError("radial_clearance must be positive and finite, got $rc"))
        isfinite(bl) && bl > 0 ||
            throw(ArgumentError("bearing_length must be positive and finite, got $bl"))
        isfinite(er) && 0 <= er < 1 ||
            throw(ArgumentError("eccentricity_ratio must be in [0, 1), got $er"))
        isfinite(α) && 0 <= α < 90 ||
            throw(ArgumentError("semicone_angle_deg must be in [0, 90)"))

        hmr = jmr + rc
        ecc = rc * er
        offset = (ecc, 0.0, 0.0)
        ming = rc - ecc
        maxg = rc + ecc

        slope = tan(deg2rad(α))
        half_bearing_length = bl * 0.5
        r0 = jmr + half_bearing_length * slope
        rL = jmr - half_bearing_length * slope
        isfinite(r0) && isfinite(rL) && min(r0, rL) > 0 ||
            throw(ArgumentError("journal end radii must be positive and finite"))
        ecc < rL ||
            throw(ArgumentError("journal eccentricity exceeds the small-end radius"))

        return new(
            jmr,
            α,
            rc,
            bl,
            er,
            hmr,
            ecc,
            offset,
            ming,
            maxg,
            slope,
            half_bearing_length,
        )
    end
end

@inline function _bearing_point_unchecked(
    p::BearingParameters,
    θ::Float64,
    z::Float64,
    ξ::Float64,
)::NTuple{3,Float64}

    Rj = p.journal_midpoint_radius + (p.half_bearing_length - z) * p.slope
    Rh = Rj + p.radial_clearance

    sinθ, cosθ = sincos(θ)
    e = p.eccentricity

    ρⱼ = e * cosθ + sqrt(Rj^2 - (e * sinθ)^2) # quadratic formula positive root
    gap = Rh - ρⱼ
    ρ = ρⱼ + (ξ * gap)

    return (ρ * cosθ, ρ * sinθ, z)
end


function bearing_point(p::BearingParameters, θ, z, ξ)
    θ, z, ξ = Float64(θ), Float64(z), Float64(ξ)

    isfinite(θ) || throw(ArgumentError("θ must be finite"))
    isfinite(z) && 0 <= z <= p.bearing_length ||
        throw(ArgumentError("z must be within the bearing length"))
    isfinite(ξ) && 0 <= ξ <= 1 || throw(ArgumentError("ξ must be in [0, 1]"))

    return _bearing_point_unchecked(p, θ, z, ξ)
end
