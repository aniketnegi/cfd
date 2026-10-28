include(joinpath(@__DIR__,"cases.jl"))
function main(output;full=false)
    fresh_output(output;suite=full ? "cylindrical-all11" : "cylindrical-representative")
    reference,_=readdlm(joinpath(@__DIR__,"reference.csv"),',';header=true)
    points=[(point[1],point[2]) for point in eachrow(reference)]
    full || filter!(point->point==(1.,.8),points)
    grids=full ? (64,128,256,512) : (64,128,256)
    mkpath(joinpath(output,"fields"))
    open(joinpath(output,"cylindrical.tsv"),"w") do io
        first=true
        for model in ("classical","production"),(aspect,ecc) in points,nt in grids
            row,field=cylindrical_case(aspect,ecc,model,nt)
            row!(io,row;header=first);first=false
            serialize(joinpath(output,"fields","$(model)-LD$(aspect)-e$(ecc)-n$(nt).jls"),field)
            @printf("%s LD=%g eps=%g nt=%d Q0=%.9f KKT=%.3g\n",model,aspect,ecc,nt,row.flow_Q0,row.complementarity)
            flush(stdout);GC.gc()
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (1,2) && (length(ARGS)==1 || ARGS[2]=="--full") || error("usage: run_cylindrical.jl FRESH_OUTPUT [--full]")
    main(ARGS[1];full="--full" in ARGS)
end
