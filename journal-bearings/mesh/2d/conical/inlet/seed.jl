# Geometry-only extraction from expts/surface-inlet-hex/mesh.jl.
# Fixed bearing/midpoint; diameter scales the local seed, not the bearing.
# No Gmsh, reservoir core, gap extrusion or pressure model.
@inline function housing_radius(p, z::Float64)
    return p.journal_midpoint_radius +
           (p.half_bearing_length - z) * p.slope +
           p.radial_clearance
end


# Convert an (x, z) point on the +Y side of the conical housing to (theta, z).
function housing_coordinates(p, x::Float64, z::Float64)
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



function reference_seed(;inlet_diameter_mm=5.)
    diameter=Float64(inlet_diameter_mm)
    isfinite(diameter) && 4<=diameter<=10 ||
        throw(ArgumentError("diameter study supports 4–10 mm only"))
    p=(journal_midpoint_radius=50., radial_clearance=.05, bearing_length=100.,
       half_bearing_length=50., slope=tan(deg2rad(10.)))
    ntheta,nz,q=128,64,8
    inlet_radius,inlet_z,patch_half_width=diameter/2,50.,1.5diameter
    curved_transition=true
    control_radius_factor,control_square_blend=2.87,.36
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


    outer_rings = Vector{Int}[rim_ring]
    for layer in 1:8
        push!(outer_rings, interpolate_ring(rim_ring, control_ring, layer / 9))
    end
    push!(outer_rings, control_ring)
    for layer in 1:9, i in 1:rim_segments
        inner,outer=outer_rings[layer],outer_rings[layer+1]
        i2=mod1(i+1,rim_segments)
        push!(surface_quads,(inner[i],outer[i],outer[i2],inner[i2]))
    end
    rim=Set(Tuple(sort([rim_ring[i],rim_ring[mod1(i+1,rim_segments)]]))
            for i in 1:rim_segments)
    R,L,ca=.05005,.1,cosd(10.)
    nodes=[(R*t,z*1e-3/ca) for (t,z) in surface_nodes]
    return assemble_surface_mesh(nodes,surface_quads,rim;R,L,ca,inlet_radius_m=inlet_radius*1e-3)
end
