# One flow per face, oriented out of its owner and into its neighbour.
# A boundary face has neighbour=0 and positive outward flow.
function face_divergence(owner, neighbour, cells)
    rows, columns, values = Int[], Int[], Float64[]
    for face in eachindex(owner)
        push!(rows,owner[face]); push!(columns,face); push!(values,1.)
        if neighbour[face] != 0
            push!(rows,neighbour[face]); push!(columns,face); push!(values,-1.)
        end
    end
    return sparse(rows,columns,values,cells,length(owner))
end

function tensor_face_operator(gt, gz, adv, ends)
    nt,nz = size(gt)
    cells = nt*nz
    count = cells+length(gz)
    owner,neighbour = zeros(Int,count),zeros(Int,count)
    boundary = fill(:internal,count)
    constant,advective = zeros(count),zeros(count)
    rows,columns,values = Int[],Int[],Float64[]
    id(i,j) = i+nt*(j-1)
    function connection(face, k, l, conductance)
        owner[face],neighbour[face] = k,l
        push!(rows,face); push!(columns,k); push!(values,conductance)
        if l != 0
            push!(rows,face); push!(columns,l); push!(values,-conductance)
        end
    end
    for j in 1:nz, i in 1:nt
        face = id(i,j)
        connection(face,id(i,j),id(mod1(i+1,nt),j),gt[i,j])
        constant[face] = advective[face] = adv[i,j]
    end
    for j in 1:nz+1, i in 1:nt
        face = cells+id(i,j)
        k = id(i,j == 1 ? 1 : j-1)
        l = j in (1,nz+1) ? 0 : id(i,j)
        connection(face,k,l,gz[i,j])
        if l == 0
            boundary[face] = j == 1 ? :large : :small
            constant[face] = -gz[i,j]*ends[j == 1 ? 1 : 2]
        end
    end
    pressure_flux = sparse(rows,columns,values,count,cells)
    divergence = face_divergence(owner,neighbour,cells)
    return (;owner,neighbour,boundary,pressure_flux,constant,advective,divergence)
end

function boundary_outflows(faces, flow)
    return (;large=sum(flow[faces.boundary .== :large]),
        small=sum(flow[faces.boundary .== :small]),
        feed=sum(flow[faces.boundary .== :feed]))
end

"""Oriented m3/s face flows and cell balances for either mesh layout.
The pressure operator satisfies q=F*p+c, A=D*F and b=-D*c, up to roundoff.
"""
function full_film_face_fluxes(a, pressure_pa)
    length(pressure_pa) == length(a.b) || throw(ArgumentError("pressure size does not match mesh"))
    all(isfinite,pressure_pa) || throw(ArgumentError("pressure must be finite"))
    face_flow = a.faces.pressure_flux*pressure_pa+a.faces.constant
    return (;face_flow,balance=a.faces.divergence*face_flow,
        boundary_outflow=boundary_outflows(a.faces,face_flow))
end
