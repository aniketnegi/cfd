# Signed Reynolds flux on the public fitted land mesh, in its metre chart.
function fitted_full_film_geometry(grid)
    all(name -> hasproperty(grid,name), (:parameters,:surface,:feed_edges,:acceptance)) ||
        throw(ArgumentError("provide a fitted grid from surface_feed_grid"))
    grid.acceptance.accepted || throw(ArgumentError("mesh has not passed acceptance"))
    p,s = grid.parameters,grid.surface
    s.geometry_mm == bearing_geometry(p) && s.R == p.housing_midpoint_radius*1e-3 &&
        s.L == p.bearing_length*1e-3 && s.ca == cosd(p.semicone_angle_deg) &&
        s.period == 2pi*s.R ||
        throw(ArgumentError("fitted chart disagrees with bearing parameters"))
    grid.quads == s.quads && grid.feed_edges == s.edge_nodes[s.boundary .== :feed] ||
        throw(ArgumentError("fitted connectivity disagrees with the surface"))
    grid.nodes == [bearing_point(p,x/s.R,
        SurfaceMeshQuality.chart_axial_position_mm(s.geometry_mm,s.ca,y),1.) for (x,y) in s.nodes] ||
        throw(ArgumentError("fitted vertices disagree with the surface"))
    n,f = length(s.quads),length(s.owner)
    size(s.centre) == (n,2) && length(s.areas) == n &&
        all(isfinite,s.centre) && all(x -> isfinite(x) && x > 0,s.areas) &&
        length(grid.areas) == length(grid.radial_gaps) == n &&
        all(x -> isfinite(x) && x > 0,grid.areas) &&
        all(x -> isfinite(x) && x > 0,grid.radial_gaps) ||
        throw(ArgumentError("invalid fitted cell geometry"))
    all(length(getproperty(s,name)) == f for name in (:neighbour,:mid,:normals,:lengths,:boundary)) &&
        all(k -> 1 <= k <= n,s.owner) && all(k -> 0 <= k <= n,s.neighbour) &&
        all(x -> isfinite(x) && x > 0,s.lengths) &&
        all(x -> all(isfinite,x),s.mid) && all(x -> all(isfinite,x) && abs(norm(x)-1) < 1e-12,s.normals) &&
        all(tag -> tag in (:internal,:large,:small,:feed),s.boundary) &&
        all(i -> (s.neighbour[i] == 0) == (s.boundary[i] in (:large,:small,:feed)),1:f) ||
        throw(ArgumentError("invalid fitted face geometry"))
    return s
end

# Weighted least-squares gradients, including the actual Dirichlet face values.
function fitted_pressure_gradients(s, boundary_values)
    n = length(s.quads)
    adjacent = [Int[] for _ in 1:n]
    for face in eachindex(s.owner)
        push!(adjacent[s.owner[face]],face)
        s.neighbour[face] == 0 || push!(adjacent[s.neighbour[face]],face)
    end
    rows,columns,gx,gy = Int[],Int[],Float64[],Float64[]
    constant = zeros(n,2)
    wrap(dx) = mod(dx+s.period/2,s.period)-s.period/2
    for cell in 1:n
        neighbours = Int[]
        displacements = Vector{Float64}[]
        for face in adjacent[cell]
            other = s.owner[face] == cell ? s.neighbour[face] : s.owner[face]
            d = other == 0 ? s.mid[face]-s.centre[cell,:] : s.centre[other,:]-s.centre[cell,:]
            d[1] = wrap(d[1])
            dot(d,d) > 0 || throw(ArgumentError("zero fitted gradient distance"))
            push!(neighbours,other); push!(displacements,d)
        end
        metric = sum(d*d'/dot(d,d) for d in displacements)
        minimum(eigvals(Symmetric(metric))) > 1e-8 || throw(ArgumentError("singular fitted gradient stencil"))
        for (face,other,d) in zip(adjacent[cell],neighbours,displacements)
            weight = metric \ (d/dot(d,d))
            push!(rows,cell); push!(columns,cell); push!(gx,-weight[1]); push!(gy,-weight[2])
            if other == 0
                constant[cell,:] .+= weight*boundary_values[face]
            else
                push!(rows,cell); push!(columns,other); push!(gx,weight[1]); push!(gy,weight[2])
            end
        end
    end
    return (;x=sparse(rows,columns,gx,n,n),y=sparse(rows,columns,gy,n,n),constant)
end

function fitted_face_operator(s, p; rpm, viscosity_pa_s, eccentricity_angle_rad, boundary_pressure,
        gap_correction_m=nothing)
    n,f = length(s.quads),length(s.owner)
    boundary_values = [s.neighbour[i] == 0 ? Float64(boundary_pressure(s.mid[i]...,s.boundary[i])) : 0. for i in 1:f]
    all(isfinite,boundary_values) || throw(ArgumentError("boundary pressure must be finite"))
    gradient = fitted_pressure_gradients(s,boundary_values)
    rows,columns,values = Int[],Int[],Float64[]
    erows,ecolumns,evalues = Int[],Int[],Float64[]
    correction_x,correction_y = zeros(f),zeros(f)
    constant,advective = zeros(f),zeros(f)
    wrap(dx) = mod(dx+s.period/2,s.period)-s.period/2
    for face in 1:f
        k,l = s.owner[face],s.neighbour[face]
        x,normal,width = s.mid[face],s.normals[face],s.lengths[face]
        z = x[2]*s.ca
        rj,rh = film_radii(p,z)
        h = corrected_normal_gap(p,x[1]/s.R,z,eccentricity_angle_rad,gap_correction_m)
        mobility = h^3/(12viscosity_pa_s)
        # In (R_ref*theta,z/cos(alpha)), diffusion is K*diag(R_ref/rh,rh/R_ref).
        metric = [s.R/rh,rh/s.R]
        d = l == 0 ? x-s.centre[k,:] : s.centre[l,:]-s.centre[k,:]
        d[1] = wrap(d[1])
        dot(normal,d) > 0 || throw(ArgumentError("nonpositive fitted pressure connection"))
        coefficient = dot(normal,metric.*normal)/dot(normal,d)
        conductance = width*mobility*coefficient
        push!(rows,face); push!(columns,k); push!(values,conductance)
        if l == 0
            constant[face] = -conductance*boundary_values[face]
        else
            push!(rows,face); push!(columns,l); push!(values,-conductance)
        end
        correction = -width*mobility*(metric.*normal-coefficient*d)
        correction_x[face],correction_y[face] = correction
        push!(erows,face); push!(ecolumns,k); push!(evalues,l == 0 ? 1. : .5)
        if l != 0
            push!(erows,face); push!(ecolumns,l); push!(evalues,.5)
        end
        advective[face] = rpm*2pi/60*rj*h/2*normal[1]*width
    end
    interpolate = sparse(erows,ecolumns,evalues,f,n)
    pressure_flux = sparse(rows,columns,values,f,n)+
        spdiagm(0=>correction_x)*interpolate*gradient.x+
        spdiagm(0=>correction_y)*interpolate*gradient.y
    constant .+= correction_x.*(interpolate*gradient.constant[:,1])+
        correction_y.*(interpolate*gradient.constant[:,2])+advective
    all(isfinite,nonzeros(pressure_flux)) && all(isfinite,constant) ||
        throw(ArgumentError("nonfinite fitted film coefficients"))
    divergence = face_divergence(s.owner,s.neighbour,n)
    return (;owner=s.owner,neighbour=s.neighbour,boundary=s.boundary,
        pressure_flux,constant,advective,divergence)
end

function assemble_fitted_full_film(grid; rpm, viscosity_pa_s, end_pressure_pa,
        eccentricity_angle_rad, feed_pressure_pa, gap_correction_m=nothing)
    isnothing(feed_pressure_pa) && throw(ArgumentError("fitted full film requires feed_pressure_pa"))
    feed = Float64(feed_pressure_pa)
    isfinite(feed) || throw(ArgumentError("feed pressure must be finite"))
    s = fitted_full_film_geometry(grid)
    p = grid.parameters
    boundary_pressure(x,y,tag) = tag == :feed ? feed : end_pressure_pa[tag == :large ? 1 : 2]
    faces = fitted_face_operator(s,p;rpm,viscosity_pa_s,eccentricity_angle_rad,boundary_pressure,gap_correction_m)
    A = faces.divergence*faces.pressure_flux
    b = -(faces.divergence*faces.constant)
    normal_gap_m = [corrected_normal_gap(p,s.centre[k,1]/s.R,s.centre[k,2]*s.ca,eccentricity_angle_rad,gap_correction_m) for k in axes(s.centre,1)]
    all(x -> isfinite(x) && x > 0,normal_gap_m) || throw(ArgumentError("invalid fitted normal gaps"))
    radial_gaps = normal_gap_m.*(1e3/s.ca)
    result_grid = eccentricity_angle_rad == get(grid,:eccentricity_angle_rad,0.) &&
        radial_gaps == grid.radial_gaps ? grid : (;grid...,radial_gaps,eccentricity_angle_rad)
    return (;A,b,grid=result_grid,normal_gap_m,faces,end_pressure_pa,feed_pressure_pa=feed,
        eccentricity_angle_rad)
end

# Integrate the curved housing surface, not planar normals of the 3D vertex quads.
function fitted_pressure_force(grid, pressure_pa)
    s = fitted_full_film_geometry(grid)
    length(pressure_pa) == length(s.quads) || throw(ArgumentError("pressure size does not match mesh"))
    all(isfinite,pressure_pa) || throw(ArgumentError("pressure must be finite"))
    p = grid.parameters
    points,weights = (-sqrt(3/5),0.,sqrt(3/5)),(5/9,8/9,5/9)
    force = zeros(3)
    for (cell,xy) in enumerate(s.localnodes), (u,wu) in zip(points,weights), (v,wv) in zip(points,weights)
        shape = [(1-u)*(1-v),(1+u)*(1-v),(1+u)*(1+v),(1-u)*(1+v)]/4
        du = [-(1-v),1-v,1+v,-(1+v)]/4
        dv = [-(1-u),-(1+u),1+u,1-u]/4
        x,xu,xv = xy'*shape,xy'*du,xy'*dv
        jacobian = xu[1]*xv[2]-xu[2]*xv[1]
        jacobian > 0 || throw(ArgumentError("inverted fitted force quadrature"))
        theta,z = x[1]/s.R,x[2]*s.ca
        _,rh = film_radii(p,z)
        weight = pressure_pa[cell]*rh/s.R*jacobian*wu*wv
        force .+= weight.*[s.ca*cos(theta),s.ca*sin(theta),sind(p.semicone_angle_deg)]
    end
    return Tuple(force)
end
