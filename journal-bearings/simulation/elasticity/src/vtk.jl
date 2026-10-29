# Reference solid coordinates in metres; displacements are a separate vector.
function write_elasticity_vtk(path,a,state)
    ispath(path) && throw(ArgumentError("refusing to overwrite $path"))
    length(state.displacement_m)==size(a.K,1) || throw(ArgumentError("solid displacement size mismatch"))
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"# vtk DataFile Version 3.0\nLinear elastic solid (SI)\nASCII\nDATASET UNSTRUCTURED_GRID")
        println(io,"POINTS $(size(a.nodes,1)) double")
        for i in axes(a.nodes,1);println(io,join(a.nodes[i,:],' '));end
        println(io,"CELLS $(length(a.cells)) $(9length(a.cells))")
        for cell in a.cells;println(io,"8 ",join(collect(cell).-1,' '));end
        println(io,"CELL_TYPES $(length(a.cells))")
        for _ in a.cells;println(io,12);end
        println(io,"POINT_DATA $(size(a.nodes,1))")
        for (name,data) in (("displacement_m",state.displacement_m),("support_reaction_n",state.reaction_n),("applied_load_n",state.load_n))
            println(io,"VECTORS $name double")
            for i in axes(a.nodes,1);println(io,join(data[3i-2:3i],' '));end
        end
        println(io,"SCALARS fixed_displacement int 1\nLOOKUP_TABLE default")
        fixed=Set(a.fixed)
        for i in axes(a.nodes,1);println(io,all(k->3*(i-1)+k in fixed,1:3) ? 1 : 0);end
    end
    return path
end
