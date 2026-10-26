using Gmsh: gmsh
using SparseArrays

include(joinpath(@__DIR__, "..", "..", "bearing.jl"))


@inline function housing_radius(p::BearingParameters, z::Float64)
    return p.journal_midpoint_radius +
           (p.half_bearing_length - z) * p.slope +
           p.radial_clearance
end


# Convert an (x, z) point on the +Y side of the conical housing to (theta, z).
function housing_coordinates(p::BearingParameters, x::Float64, z::Float64)
    radius = housing_radius(p, z)
    abs(x) < radius || throw(ArgumentError("inlet point lies outside housing"))
    y = sqrt(radius^2 - x^2)
    return (atan(y, x), z)
end


# Return the 4q unique vertices around a q-by-q rectangle, counterclockwise.
function rectangle_perimeter(ids::Matrix{Int}, i0::Int, j0::Int, q::Int)
    result = Int[]
    append!(result, (ids[i0 + t, j0] for t in 0:(q - 1)))
    append!(result, (ids[i0 + q, j0 + t] for t in 0:(q - 1)))
    append!(result, (ids[i0 + q - t, j0 + q] for t in 0:(q - 1)))
    append!(result, (ids[i0, j0 + q - t] for t in 0:(q - 1)))
    return result
end


"""Place a 2q-cell support at a fixed physical width, with linear spacing outside."""
function patch_axis(count, centre, extent, position, half_width, q)
    lo, hi = centre - q, centre + q
    0 < lo < hi < count || throw(ArgumentError("support must lie inside grid"))
    0 < position - 2half_width < position + 2half_width < extent ||
        throw(ArgumentError("physical support crosses domain boundary"))
    axis = [i <= lo ? (position - 2half_width) * i / lo :
            i <= hi ? position + 2half_width * (i - centre) / q :
            position + 2half_width + (extent - position - 2half_width) *
                (i - hi) / (count - hi) for i in 0:count]
    axis[1], axis[end] = 0.0, extent
    return axis
end


"""Solve the four-neighbour Laplace stencil between fixed q and 2q perimeters."""
function harmonic_support!(nodes, grid_ids, p, ic, jc, q)
    unknown = [(i, j) for j in (jc-q+1):(jc+q-1), i in (ic-q+1):(ic+q-1)
        if !(abs(i-ic) <= q÷2 && abs(j-jc) <= q÷2)]
    row_for = Dict(index => row for (row, index) in enumerate(unknown))
    rows, cols, values = Int[], Int[], Float64[]
    rhs = zeros(length(unknown), 2)
    for (row, (i, j)) in enumerate(unknown)
        push!(rows, row); push!(cols, row); push!(values, 4.0)
        for index in ((i-1, j), (i+1, j), (i, j-1), (i, j+1))
            other = get(row_for, index, 0)
            if other == 0
                theta, z = nodes[grid_ids[index...]]
                rhs[row, 1] += housing_radius(p, z) * cos(theta)
                rhs[row, 2] += z
            else
                push!(rows, row); push!(cols, other); push!(values, -1.0)
            end
        end
    end
    matrix = sparse(rows, cols, values, length(unknown), length(unknown))
    solution = matrix \ rhs
    for (row, index) in enumerate(unknown)
        nodes[grid_ids[index...]] = housing_coordinates(p, solution[row, 1], solution[row, 2])
    end
    return nothing
end


"""
    bearing_surface(p; ntheta, nz, inlet_radius, inlet_z, q,
                    inner_layers, outer_layers, patch_half_width=nothing,
                    curved_transition=false, core_fraction=...,
                    control_radius_factor=2.87, control_square_blend=0.36)

Build a quadrilateral master mesh in `(theta, z)` on the housing surface.
Lengths are millimetres; `inlet_radius` is a radius, not a diameter.
The drilling cylinder is directed along +Y.
Its circular rim is filled by a square core and annular quadrilaterals, then
connected to the ordinary periodic bearing grid by transition rings.

Default construction is unchanged. `patch_half_width` fixes the replaced patch's
axial half-width and circumferential arc half-width at the inlet centre, while
counts control resolution. Its outer support has twice that half-width.
With `nothing`, patch dimensions still follow the uniform background grid.

`curved_transition=true` uses tangent-spaced core nodes and a circle/square
control ring, with a harmonic transition to the fixed 2q support perimeter.
`core_fraction` is core corner radius / inlet radius (default 0.9 for the curved
candidate, 0.5 otherwise). A value of 1 collapses inner-ring corner edges.
The control circle radius is `control_radius_factor * inlet_radius`; blend 0
is circular and 1 retains the background rectangle. These controls do not
guarantee positive Jacobians or a quality threshold; measure the swept Hex8 mesh.
"""
function bearing_surface(
    p::BearingParameters;
    ntheta::Int,
    nz::Int,
    inlet_radius::Float64,
    inlet_z::Float64,
    q::Int,
    inner_layers::Int,
    outer_layers::Int,
    patch_half_width::Union{Nothing,Float64} = nothing,
    curved_transition::Bool = false,
    core_fraction::Float64 = curved_transition ? 0.9 : 0.5,
    control_radius_factor::Float64 = 2.87,
    control_square_blend::Float64 = 0.36,
)
    ntheta >= 4 && ntheta % 4 == 0 ||
        throw(ArgumentError("ntheta must be divisible by four"))
    nz >= 2 && nz % 2 == 0 || throw(ArgumentError("nz must be even"))
    q >= 2 && q % 2 == 0 || throw(ArgumentError("q must be positive and even"))
    inner_layers >= 1 && outer_layers >= 1 ||
        throw(ArgumentError("O-grid layer counts must be positive"))
    isfinite(inlet_radius) && inlet_radius > 0 ||
        throw(ArgumentError("inlet_radius must be positive and finite"))
    0 < core_fraction < 1 || throw(ArgumentError("core_fraction must lie in (0, 1)"))
    isfinite(control_radius_factor) && control_radius_factor > 1 ||
        throw(ArgumentError("control_radius_factor must be finite and greater than 1"))
    0 <= control_square_blend <= 1 ||
        throw(ArgumentError("control_square_blend must lie in [0, 1]"))
    inlet_radius < inlet_z < p.bearing_length - inlet_radius ||
        throw(ArgumentError("inlet must lie between axial ends"))

    # Put inlet centre at theta=pi/2 (+Y). The replaced background patch has
    # q cells in each direction and shares its perimeter with O-grid.
    ic = ntheta ÷ 4 + 1
    axial_index = inlet_z * nz / p.bearing_length
    jc = round(Int, axial_index) + 1
    isapprox(axial_index, jc - 1; atol = 1e-12, rtol = 0) ||
        throw(ArgumentError("inlet centre must coincide with axial grid node"))
    i0 = ic - q ÷ 2
    j0 = jc - q ÷ 2
    1 <= i0 && i0 + q <= ntheta || throw(ArgumentError("inlet patch crosses seam"))
    1 <= j0 && j0 + q <= nz + 1 || throw(ArgumentError("inlet patch crosses end"))
    if curved_transition || patch_half_width !== nothing
        1 < ic-q && ic+q < ntheta || throw(ArgumentError("support crosses seam"))
        1 < jc-q && jc+q < nz+1 || throw(ArgumentError("support crosses end"))
    end

    theta_axis = [2pi * (i - 1) / ntheta for i in 1:ntheta]
    z_axis = [p.bearing_length * (j - 1) / nz for j in 1:(nz+1)]
    if patch_half_width !== nothing
        isfinite(patch_half_width) && patch_half_width > inlet_radius ||
            throw(ArgumentError("patch_half_width must be finite and exceed inlet radius"))
        theta_axis = patch_axis(ntheta, ic-1, 2pi, pi/2,
            patch_half_width / housing_radius(p, inlet_z), q)[1:end-1]
        z_axis = patch_axis(nz, jc-1, p.bearing_length, inlet_z, patch_half_width, q)
    end
    # The complete circular mouth must fit inside the original patch, including
    # the smaller-radius axial side of a conical housing.
    minimum(housing_radius(p, z) * sin(pi/2 - theta_axis[i0])
        for z in (z_axis[j0], z_axis[j0+q])) > inlet_radius &&
        min(inlet_z-z_axis[j0], z_axis[j0+q]-inlet_z) > inlet_radius ||
        throw(ArgumentError("background patch does not enclose inlet"))

    surface_nodes = NTuple{2,Float64}[]
    grid_ids = zeros(Int, ntheta, nz + 1)

    # Keep regular background nodes except unused nodes inside replaced patch.
    for j in 1:(nz + 1), i in 1:ntheta
        if !(i0 < i < i0 + q && j0 < j < j0 + q)
            theta = theta_axis[i]
            z = z_axis[j]
            push!(surface_nodes, (theta, z))
            grid_ids[i, j] = length(surface_nodes)
        end
    end

    surface_quads = NTuple{4,Int}[]
    inlet_quads = Bool[]

    # Periodic background grid, excluding q-by-q opening replaced by O-grid.
    for j in 1:nz, i in 1:ntheta
        in_patch = i0 <= i < i0 + q && j0 <= j < j0 + q
        in_patch && continue
        i2 = i == ntheta ? 1 : i + 1
        push!(surface_quads, (
            grid_ids[i, j],
            grid_ids[i2, j],
            grid_ids[i2, j + 1],
            grid_ids[i, j + 1],
        ))
        push!(inlet_quads, false)
    end

    control_ring = rectangle_perimeter(grid_ids, i0, j0, q)
    rim_segments = 4q

    if curved_transition
        for (m, id) in enumerate(control_ring)
            angle = 5pi/4 + 2pi * (m-1) / rim_segments
            theta, z = surface_nodes[id]
            x = housing_radius(p, z) * cos(theta)
            radius = control_radius_factor * inlet_radius
            x = (1-control_square_blend) * (-radius*cos(angle)) + control_square_blend*x
            z = (1-control_square_blend) * (inlet_z+radius*sin(angle)) + control_square_blend*z
            hypot(x, z-inlet_z) > inlet_radius ||
                throw(ArgumentError("control ring must lie outside inlet"))
            surface_nodes[id] = housing_coordinates(p, x, z)
        end
        harmonic_support!(surface_nodes, grid_ids, p, ic, jc, q)
    end

    half_side = core_fraction * inlet_radius / sqrt(2.0)
    central_ids = Matrix{Int}(undef, q + 1, q + 1)
    for j in 1:(q + 1), i in 1:(q + 1)
        u, v = 2(i-1)/q - 1, 2(j-1)/q - 1
        s = curved_transition ? half_side*tan(pi*u/4) : -half_side + 2half_side*(i-1)/q
        z = curved_transition ? inlet_z + half_side*tan(pi*v/4) :
            inlet_z - half_side + 2half_side*(j-1)/q
        push!(surface_nodes, housing_coordinates(p, -s, z))
        central_ids[i, j] = length(surface_nodes)
    end

    for j in 1:q, i in 1:q
        push!(surface_quads, (
            central_ids[i, j],
            central_ids[i + 1, j],
            central_ids[i + 1, j + 1],
            central_ids[i, j + 1],
        ))
        push!(inlet_quads, true)
    end
    central_ring = rectangle_perimeter(central_ids, 1, 1, q)

    # Exact intersection of conical housing and drilling cylinder
    # x^2 + (z-inlet_z)^2 = inlet_radius^2, on +Y housing side.
    rim_ring = Int[]
    for m in 0:(rim_segments - 1)
        angle = 5pi / 4 + 2pi * m / rim_segments
        s = inlet_radius * cos(angle)
        z = inlet_z + inlet_radius * sin(angle)
        push!(surface_nodes, housing_coordinates(p, -s, z))
        push!(rim_ring, length(surface_nodes))
    end

    function interpolate_ring(inner::Vector{Int}, outer::Vector{Int}, fraction::Float64)
        ring = Int[]
        for (inner_id, outer_id) in zip(inner, outer)
            theta1, z1 = surface_nodes[inner_id]
            theta2, z2 = surface_nodes[outer_id]
            x1 = housing_radius(p, z1) * cos(theta1)
            x2 = housing_radius(p, z2) * cos(theta2)
            x = (1 - fraction) * x1 + fraction * x2
            z = (1 - fraction) * z1 + fraction * z2
            push!(surface_nodes, housing_coordinates(p, x, z))
            push!(ring, length(surface_nodes))
        end
        return ring
    end

    inner_rings = Vector{Int}[central_ring]
    for layer in 1:(inner_layers - 1)
        push!(inner_rings, interpolate_ring(central_ring, rim_ring, layer / inner_layers))
    end
    push!(inner_rings, rim_ring)

    outer_rings = Vector{Int}[rim_ring]
    for layer in 1:(outer_layers - 1)
        push!(outer_rings, interpolate_ring(rim_ring, control_ring, layer / outer_layers))
    end
    push!(outer_rings, control_ring)

    function append_ring_quads!(rings::Vector{Vector{Int}}, is_inlet::Bool)
        for layer in 1:(length(rings) - 1)
            inner = rings[layer]
            outer = rings[layer + 1]
            for i in 1:rim_segments
                i2 = i == rim_segments ? 1 : i + 1
                push!(surface_quads, (inner[i], outer[i], outer[i2], inner[i2]))
                push!(inlet_quads, is_inlet)
            end
        end
    end

    append_ring_quads!(inner_rings, true)
    append_ring_quads!(outer_rings, false)

    return surface_nodes, surface_quads, inlet_quads
end


"""Map every housing-surface vertex through all clearance layers."""
function bearing_nodes(
    p::BearingParameters,
    surface_nodes::Vector{NTuple{2,Float64}},
    nxi::Int,
)
    nxi >= 1 || throw(ArgumentError("nxi must be positive"))
    points = Matrix{NTuple{3,Float64}}(undef, length(surface_nodes), nxi + 1)
    for k in 1:(nxi + 1)
        xi = k == nxi + 1 ? 1.0 : (k - 1) / nxi
        for v in eachindex(surface_nodes)
            theta, z = surface_nodes[v]
            points[v, k] = bearing_point(p, theta, z, xi)
        end
    end
    return points
end


"""Sweep each master-surface Quad4 through clearance into `nxi` Hex8 cells."""
function bearing_hexes(
    surface_quads::Vector{NTuple{4,Int}},
    node_count::Int,
    nxi::Int,
)
    ids = LinearIndices((node_count, nxi + 1))
    hexes = Matrix{NTuple{8,Int}}(undef, length(surface_quads), nxi)
    for k in 1:nxi, q in eachindex(surface_quads)
        a, b, c, d = surface_quads[q]
        hexes[q, k] = (
            ids[a, k], ids[b, k], ids[c, k], ids[d, k],
            ids[a, k + 1], ids[b, k + 1], ids[c, k + 1], ids[d, k + 1],
        )
    end
    return hexes
end


"""Extract outward Quad4 faces for five external boundary groups."""
function bearing_quads(
    p::BearingParameters,
    surface_nodes::Vector{NTuple{2,Float64}},
    surface_quads::Vector{NTuple{4,Int}},
    inlet_quads::Vector{Bool},
    hexes::Matrix{NTuple{8,Int}},
)
    journal = NTuple{4,Int}[]
    housing = NTuple{4,Int}[]
    inlet = NTuple{4,Int}[]
    z0 = NTuple{4,Int}[]
    zL = NTuple{4,Int}[]
    nxi = size(hexes, 2)
    side_corners = ((1, 2, 6, 5), (2, 3, 7, 6),
                    (3, 4, 8, 7), (4, 1, 5, 8))

    for q in eachindex(surface_quads)
        first_hex = hexes[q, 1]
        push!(journal, (first_hex[1], first_hex[4], first_hex[3], first_hex[2]))

        last_hex = hexes[q, nxi]
        outer = (last_hex[5], last_hex[6], last_hex[7], last_hex[8])
        push!(inlet_quads[q] ? inlet : housing, outer)

        a, b, c, d = surface_quads[q]
        surface_sides = ((a, b), (b, c), (c, d), (d, a))
        for side in 1:4
            u, v = surface_sides[side]
            zu = surface_nodes[u][2]
            zv = surface_nodes[v][2]
            output = zu == zv == 0.0 ? z0 :
                     zu == zv == p.bearing_length ? zL : nothing
            output === nothing && continue
            corners = side_corners[side]
            for k in 1:nxi
                h = hexes[q, k]
                push!(output, ntuple(i -> h[corners[i]], 4))
            end
        end
    end

    return (; journal, housing, inlet, z0, zL)
end


function write_bearing_mesh(path, points, hexes, boundaries)
    path = abspath(path)
    mkpath(dirname(path))
    gmsh.initialize(String[], false)
    try
        gmsh.model.add("surface-inlet-bearing")
        for entity in 2:6
            gmsh.model.addDiscreteEntity(2, entity)
        end
        gmsh.model.addDiscreteEntity(3, 1, collect(2:6))
        gmsh.model.mesh.addNodes(
            3, 1, collect(1:length(points)), collect(Iterators.flatten(points)),
        )
        gmsh.model.mesh.addElementsByType(
            1, 5, collect(1:length(hexes)), collect(Iterators.flatten(hexes)),
        )
        gmsh.model.addPhysicalGroup(3, [1], 1, "fluid")

        next_element = length(hexes) + 1
        groups = ((2, "journal", boundaries.journal),
                  (3, "housing", boundaries.housing),
                  (4, "inlet", boundaries.inlet),
                  (5, "z0", boundaries.z0),
                  (6, "zL", boundaries.zL))
        for (entity, name, quads) in groups
            tags = collect(next_element:(next_element + length(quads) - 1))
            gmsh.model.mesh.addElementsByType(
                entity, 3, tags, collect(Iterators.flatten(quads)),
            )
            gmsh.model.addPhysicalGroup(2, [entity], entity, name)
            next_element += length(quads)
        end

        gmsh.model.mesh.reclassifyNodes()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.option.setNumber("Mesh.PreserveNumberingMsh2", 1)
        gmsh.write(path)
    finally
        gmsh.finalize()
    end
    return path
end


if abspath(PROGRAM_FILE) == @__FILE__
    # One review mesh. Lengths are millimetres.
    p = BearingParameters(
        journal_midpoint_radius = 50.0,
        semicone_angle_deg = 10.0,
        radial_clearance = 0.05,
        bearing_length = 100.0,
        eccentricity_ratio = 0.6,
    )
    inlet_diameter = 0.1 * p.journal_midpoint_radius

    surface_nodes, surface_quads, inlet_quads = bearing_surface(
        p;
        ntheta = 256,
        nz = 100,
        inlet_radius = inlet_diameter / 2,
        inlet_z = 50.0,
        q = 8,
        inner_layers = 2,
        outer_layers = 4,
        patch_half_width = 3 * inlet_diameter / 2,
    )
    points = bearing_nodes(p, surface_nodes, 4)
    hexes = bearing_hexes(surface_quads, length(surface_nodes), 4)
    boundaries = bearing_quads(p, surface_nodes, surface_quads, inlet_quads, hexes)
    path = joinpath(dirname(@__DIR__), "out", "surface-inlet-hex", "bearing.msh")
    write_bearing_mesh(path, points, hexes, boundaries)
    println("Mesh: $path")
    println("$(length(points)) nodes, $(length(hexes)) hexes, " *
            "$(length(boundaries.inlet)) inlet faces")
end
