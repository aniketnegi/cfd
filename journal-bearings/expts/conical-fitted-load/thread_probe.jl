# Post-run diagnostic only: no mesh generation or pressure solve.
include(joinpath(@__DIR__,"solve.jl"))
using Libdl
BLAS.set_num_threads(1)
println("Julia threads: ",Threads.nthreads())
println("Julia BLAS interface: ",BLAS.get_num_threads()," ",BLAS.get_config())
for path in filter(p->occursin("openblas",lowercase(p)),Libdl.dllist())
    println("Loaded: ",path)
    handle=Libdl.dlopen(path)
    for name in (:openblas_get_num_threads,:openblas_get_num_threads64_)
        pointer=Libdl.dlsym_e(handle,name)
        pointer==C_NULL || println("  ",name," => ",ccall(pointer,Cint,()))
    end
end
