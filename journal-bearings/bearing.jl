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

"""
    bearing_geometry(p)

Pose-independent lengths in mm and the dimensionless cached slope for surface
measurements. Keep the half-length and slope rather than deriving them again
after unit conversion.
"""
function bearing_geometry(p::BearingParameters)
    return (;
        journal_midpoint_radius=p.journal_midpoint_radius,
        bearing_length=p.bearing_length,
        half_bearing_length=p.half_bearing_length,
        radial_clearance=p.radial_clearance,
        slope=p.slope,
    )
end

"""
    journal_radius(geometry, z)
    housing_radius(geometry, z)

Analytic cone radii in mm, about their respective axes. Accept BearingParameters
or its bearing_geometry record. These low-level functions extrapolate in z;
the point and radial-geometry interfaces below check the bearing bounds.
"""
@inline journal_radius(geometry, z::Real) =
    geometry.journal_midpoint_radius + (geometry.half_bearing_length - z) * geometry.slope

@inline housing_radius(geometry, z::Real) =
    journal_radius(geometry, z) + geometry.radial_clearance

@inline function _bearing_radial_geometry(p::BearingParameters, θ::Float64, z::Float64)
    radius = journal_radius(p, z)
    housing = housing_radius(p, z)
    sin_theta, cos_theta = sincos(θ)
    journal_ray = p.eccentricity * cos_theta +
        sqrt(radius^2 - (p.eccentricity * sin_theta)^2)
    return (;
        journal_ray_radius=journal_ray,
        housing_radius=housing,
        radial_gap=housing-journal_ray,
        sin_theta,
        cos_theta,
    )
end

function _checked_bearing_coordinates(p::BearingParameters, θ, z)
    θ, z = Float64(θ), Float64(z)
    isfinite(θ) || throw(ArgumentError("θ must be finite"))
    isfinite(z) && 0 <= z <= p.bearing_length ||
        throw(ArgumentError("z must be within the bearing length"))
    return θ, z
end

"""
    bearing_radial_geometry(p, theta, z)

Evaluate the journal ray intersection, housing radius and radial gap in mm.
Theta is in radians and z must lie within the bearing. This is the same radial
geometry used by bearing_point and both surface mesh constructors.
"""
function bearing_radial_geometry(p::BearingParameters, θ, z)
    θ, z = _checked_bearing_coordinates(p, θ, z)
    return _bearing_radial_geometry(p, θ, z)
end

@inline function _bearing_point_unchecked(
    p::BearingParameters,
    θ::Float64,
    z::Float64,
    ξ::Float64,
)::NTuple{3,Float64}

    # The housing is fixed: do not recover it by cancelling the journal
    # displacement through the gap, or evaluate an unused ray intersection.
    if ξ == 1
        sin_theta, cos_theta = sincos(θ)
        radius = housing_radius(p, z)
        return (radius * cos_theta, radius * sin_theta, z)
    end

    geometry = _bearing_radial_geometry(p, θ, z)
    radius = geometry.journal_ray_radius + ξ * geometry.radial_gap
    return (radius * geometry.cos_theta, radius * geometry.sin_theta, z)
end


function bearing_point(p::BearingParameters, θ, z, ξ)
    θ, z = _checked_bearing_coordinates(p, θ, z)
    ξ = Float64(ξ)
    isfinite(ξ) && 0 <= ξ <= 1 || throw(ArgumentError("ξ must be in [0, 1]"))

    return _bearing_point_unchecked(p, θ, z, ξ)
end
