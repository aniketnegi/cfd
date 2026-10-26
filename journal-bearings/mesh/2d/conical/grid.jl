# Shared periodic circumferential–axial Quad4 grid on the housing cone.
# Lengths are mm, theta is in radians. This grid defines no feed or solver conditions.
# Caller loads bearing.jl; construction does not load exporters or solvers.

function housing_grid(p::BearingParameters, n_theta::Int, n_axial::Int)
    n_theta >= 3 || throw(ArgumentError("n_theta must be at least 3"))
    n_axial >= 1 || throw(ArgumentError("n_axial must be positive"))
    p.eccentricity < p.journal_midpoint_radius - p.half_bearing_length * p.slope ||
        throw(ArgumentError("journal eccentricity exceeds the small-end radius"))

    dtheta = 2pi / n_theta
    dz = p.bearing_length / n_axial
    ds = dz / cosd(p.semicone_angle_deg)
    node_id(i, j) = i + n_theta * (j - 1)
    nodes = Vector{NTuple{3,Float64}}(undef, n_theta * (n_axial + 1))
    quads = Vector{NTuple{4,Int}}(undef, n_theta * n_axial)
    radial_gaps = Vector{Float64}(undef, length(quads))
    areas = Vector{Float64}(undef, length(quads))

    # z=0 is the large end. Unique circumferential nodes wrap at the seam.
    for j in 1:(n_axial + 1)
        z = (j - 1) * dz
        for i in 1:n_theta
            theta = (i - 1) * dtheta
            nodes[node_id(i, j)] = bearing_point(p, theta, z, 1.0)
        end
    end

    for j in 1:n_axial
        z = (j - 0.5) * dz
        for i in 1:n_theta
            next_i = i == n_theta ? 1 : i + 1
            k = i + n_theta * (j - 1)
            quads[k] = (node_id(i, j), node_id(next_i, j),
                        node_id(next_i, j + 1), node_id(i, j + 1))

            theta = (i - 0.5) * dtheta
            xh, yh, _ = bearing_point(p, theta, z, 1.0)
            xj, yj, _ = bearing_point(p, theta, z, 0.0)
            radial_gaps[k] = hypot(xh, yh) - hypot(xj, yj)
            # Analytic cone area, not the area of the faceted output quad.
            areas[k] = hypot(xh, yh) * dtheta * ds
        end
    end

    minimum(radial_gaps) > 0 || error("journal intersects housing")
    return (; nodes, quads, radial_gaps, areas, n_theta, n_axial)
end
