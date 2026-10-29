# Fixed nominal position, isothermal JFO and small-strain elastic solids.
# No load-control loop, contact, inertia or moving-mesh kinematics.
isdefined(@__MODULE__,:solve_conical_jfo) || include(joinpath(@__DIR__,"../../reynolds/src/conical_jfo.jl"))
isdefined(@__MODULE__,:assemble_housing) || include(joinpath(@__DIR__,"../../elasticity/src/housing.jl"))
isdefined(@__MODULE__,:assemble_journal) || include(joinpath(@__DIR__,"../../elasticity/src/journal.jl"))
include(joinpath(@__DIR__,"interface.jl"))
fsi_map(f,nt)=NamedTuple{keys(nt)}(Tuple(f(v) for v in values(nt)))

function prepare_housing_fsi(grid;youngs_modulus_pa,poisson_ratio,wall_thickness_m,
        n_radial,feed_diameter_m=nothing,feed_z_m=nothing)
    solid=assemble_housing(grid;youngs_modulus_pa,poisson_ratio,wall_thickness_m,
        n_radial,feed_diameter_m,feed_z_m)
    gap=housing_gap_operator(solid)
    housing=(;assembly=solid,gap_matrix=gap.matrix,opening_sign=1,
        interface_nodes=1:solid.mesh.base_nodes,reservoir_load_per_pa=solid.passage_load_per_pa)
    return (;solid,gap,grid,bodies=(;housing))
end
function add_elastic_journal(model;youngs_modulus_pa,poisson_ratio,n_theta,n_axial,n_radial)
    haskey(model.bodies,:journal) && throw(ArgumentError("model already has an elastic journal"))
    a=assemble_journal(model.grid;youngs_modulus_pa,poisson_ratio,n_theta,n_axial,n_radial)
    journal=(;assembly=a,gap_matrix=journal_gap_operator(a,model.gap.points),opening_sign=-1,
        interface_nodes=a.mesh.outer_nodes,reservoir_load_per_pa=a.reservoir_load_per_pa)
    return (;model...,bodies=(;model.bodies...,journal))
end
function fsi_pressure_reference(model,fluid_options,external_pressure_pa)
    isfinite(external_pressure_pa) || throw(ArgumentError("external pressure must be finite"))
    haskey(fluid_options,:gap_correction_m) && throw(ArgumentError("FSI owns gap_correction_m"))
    haskey(fluid_options,:feed_cells) && throw(ArgumentError("FSI uses end supply or an actual rim, not reservoir cells"))
    return hasproperty(model.grid,:feed_edges) ? fluid_options.feed_pressure_pa-external_pressure_pa : 0.
end
fsi_assembly_options(options)=(;(k=>v for (k,v) in pairs(options) if k ∉ (:max_iterations,:tolerance))...)
function fsi_gap_values(model,u)
    parts=Tuple(begin
        body=model.bodies[key];value=body.gap_matrix*u[key]
        body.opening_sign>0 ? value : -value
    end for key in keys(model.bodies))
    return reduce(+,parts)
end
fsi_gap_field(model,u)=HousingGapField(model.gap.index,fsi_gap_values(model,u))
function fsi_gap_lower_bound(model,u::NamedTuple)
    motion=sum(maximum(norm(view(u[key],3i-2:3i)) for i in model.bodies[key].interface_nodes) for key in keys(model.bodies))
    p=model.grid.parameters
    return p.minimum_gap*1e-3*cosd(p.semicone_angle_deg)-motion
end
function fsi_elastic_targets(model,pressure,passage)
    return fsi_map(model.bodies) do body
        a=body.assembly;land_load=a.pressure_operator*pressure
        load=land_load+body.reservoir_load_per_pa*passage
        state=solve_elasticity(a,load)
        (;state...,load_n=load,land_load_n=land_load,
            normal_displacement_cell_mean_m=(a.pressure_operator'*state.displacement_m)./a.interface_areas_m2)
    end
end
function fsi_stresses(model,states)
    return NamedTuple{keys(model.bodies)}(Tuple(elastic_stress(model.bodies[k].assembly,states[k].displacement_m) for k in keys(model.bodies)))
end
function solve_one_way_fsi(model;fluid_options::NamedTuple,external_pressure_pa::Real=0.)
    passage=fsi_pressure_reference(model,fluid_options,external_pressure_pa)
    fluid=solve_conical_jfo(model.grid;fluid_options...)
    solids=fsi_elastic_targets(model,fluid.pressure_pa.-external_pressure_pa,passage)
    u=fsi_map(state->state.displacement_m,solids)
    correction=fsi_gap_field(model,u)
    deformed=assemble_conical_jfo(model.grid;fsi_assembly_options(fluid_options)...,gap_correction_m=correction)
    stresses=fsi_stresses(model,solids)
    return (;mode=:one_way,fluid,solids,stresses,solid=solids.housing,stress=stresses.housing,
        gap_correction_at_assembly_points_m=correction.values,
        predicted_deformed_normal_gap_m=deformed.normal_gap_m,
        minimum_predicted_gap_m=minimum(deformed.normal_gap_m),
        conservative_gap_lower_bound_m=fsi_gap_lower_bound(model,u),
        pressure_recomputed_on_deformation=false)
end

function solve_two_way_fsi(model;fluid_options::NamedTuple,external_pressure_pa::Real=0.,
        relaxation::Real=.5,max_iterations::Int=100,gap_tolerance_m::Real=1e-10,
        pressure_tolerance::Real=1e-5,solid_tolerance::Real=1e-5,observer=nothing)
    isfinite(relaxation) && 0<relaxation<=1 || throw(ArgumentError("invalid coupling relaxation"))
    max_iterations>0 || throw(ArgumentError("invalid coupling iteration cap"))
    all(x->isfinite(x)&&x>0,(gap_tolerance_m,pressure_tolerance,solid_tolerance)) || throw(ArgumentError("invalid coupling tolerance"))
    passage=fsi_pressure_reference(model,fluid_options,external_pressure_pa)
    u=fsi_map(body->zeros(size(body.assembly.K,1)),model.bodies)
    previous_pressure=nothing;history=NamedTuple[]
    for iteration in 1:max_iterations
        correction=fsi_gap_field(model,u)
        fluid=solve_conical_jfo(model.grid;fluid_options...,gap_correction_m=correction)
        targets=fsi_elastic_targets(model,fluid.pressure_pa.-external_pressure_pa,passage)
        target_u=fsi_map(state->state.displacement_m,targets)
        du=NamedTuple{keys(u)}(Tuple(target_u[k]-u[k] for k in keys(u)))
        defect=norm(fsi_gap_values(model,du),Inf)
        residuals=NamedTuple{keys(u)}(Tuple(model.bodies[k].assembly.K*u[k]-targets[k].load_n for k in keys(u)))
        relative=NamedTuple{keys(u)}(Tuple(norm(residuals[k][model.bodies[k].assembly.free],Inf)/max(norm(targets[k].load_n,Inf),1e-30) for k in keys(u)))
        relative_solid=maximum(values(relative))
        change=isnothing(previous_pressure) ? Inf : norm(fluid.pressure_above_cavitation_pa-previous_pressure)/max(norm(fluid.pressure_above_cavitation_pa),1e-30)
        row=(;iteration,gap_defect_m=defect,pressure_relative_change=change,
            solid_relative_residual=relative_solid,minimum_gap_m=minimum(fluid.normal_gap_m),
            peak_pressure_pa=maximum(fluid.pressure_pa),minimum_fill=minimum(fluid.fill_fraction),
            liquid_global_balance=fluid.diagnostics.relative_global_imbalance)
        push!(history,row)
        isnothing(observer) || observer(row)
        if defect<=gap_tolerance_m && change<=pressure_tolerance && relative_solid<=solid_tolerance
            solids=NamedTuple{keys(u)}(Tuple(begin
                a=model.bodies[k].assembly;reaction=zeros(length(u[k]));reaction[a.fixed]=residuals[k][a.fixed]
                (;displacement_m=u[k],reaction_n=reaction,load_n=targets[k].load_n,
                    land_load_n=targets[k].land_load_n,relative_residual=relative[k],
                    strain_energy_j=dot(u[k],a.K*u[k])/2,
                    normal_displacement_cell_mean_m=(a.pressure_operator'*u[k])./a.interface_areas_m2)
            end for k in keys(u)))
            stresses=fsi_stresses(model,solids)
            return (;mode=:two_way,fluid,solids,stresses,solid=solids.housing,stress=stresses.housing,history,
                gap_correction_at_assembly_points_m=correction.values,
                conservative_gap_lower_bound_m=fsi_gap_lower_bound(model,u),
                pressure_recomputed_on_deformation=true)
        end
        previous_pressure=fluid.pressure_above_cavitation_pa
        u=NamedTuple{keys(u)}(Tuple((1-relaxation).*u[k].+relaxation.*target_u[k] for k in keys(u)))
        yield()
    end
    error("two-way FSI did not converge after $max_iterations iterations; last state=$(last(history))")
end
