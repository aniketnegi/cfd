include(joinpath(@__DIR__,"spectral.jl"))
using TOML, SHA
function profiles(directory)
    ispath(directory) && error("refusing overwrite")
    mkpath(directory)
    files=[joinpath(@__DIR__,name) for name in
        ("reference_profiles.jl","spectral.jl","table-6-5.csv","reference.toml")]
    open(joinpath(directory,"inputs.toml"),"w") do io
        TOML.print(io,Dict("suite"=>"independent-spectral-profiles","julia_version"=>string(VERSION),
            "modes"=>320,"L_D"=>1.,"physical_validation"=>false,
            "source_sha256"=>Dict(basename(path)=>bytes2hex(sha256(read(path))) for path in files));sorted=true)
    end
    _,header=readdlm(joinpath(@__DIR__,"table-6-5.csv"),',';header=true)
    eccentricities=parse.(Float64,String.(vec(header)[2:end]))
    theta=collect(range(0.,2pi;length=8193))
    angular=sin.(theta .* transpose(1:320))
    w=collect(range(-1.,1.;length=1001))
    open(joinpath(directory,"spectral-peaks.tsv"),"w") do io
        println(io,"epsilon\tmaximum_P_midplane\ttheta_peak_deg")
        for ecc in eccentricities
            basis=spectral_basis(ecc,320)
            coefficients=basis.vectors*(basis.amplitudes.*[1-cosh_ratio(a,0.,1.) for a in basis.alpha])
            pressure=angular*coefficients
            index=argmax(pressure)
            writedlm(joinpath(directory,"LD1-eps$(ecc)-circumferential.tsv"),
                vcat(["theta_deg" "P_midplane"],hcat(rad2deg.(theta),pressure)),'\t')
            at160=reference_field_at_angle(basis,deg2rad(160.),w)
            writedlm(joinpath(directory,"LD1-eps$(ecc)-axial.tsv"),
                vcat(["w_from_midplane" "P_at_160_deg"],hcat(w,at160)),'\t')
            println(io,join((ecc,pressure[index],rad2deg(theta[index])),'\t'))
        end
    end
end
function reference_field_at_angle(basis,theta,w)
    angular=transpose(sin.(collect(1:320).*theta))*basis.vectors
    axial=[basis.amplitudes[n]*(1-cosh_ratio(basis.alpha[n],z,1.)) for n in 1:320,z in w]
    return vec(angular*axial)
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("usage: reference_profiles.jl FRESH_OUTPUT")
    profiles(ARGS[1])
end
