# Bilinear point sampling of a periodic tensor field at physical coordinates.
function sample_pressure(theta,z,pressure,query_theta,query_z)
    nt,nz=length(theta),length(z)
    matrix=reshape(pressure,nt,nz)
    t=mod(query_theta-first(theta),2pi)+first(theta)
    i=searchsortedlast(theta,t);next=mod1(i+1,nt)
    high=next==1 ? first(theta)+2pi : theta[next]
    alpha=(t-theta[i])/(high-theta[i])
    j=searchsortedlast(z,query_z)
    1<=j<=nz || error("axial sample outside pressure-centre span")
    j==nz && query_z!=z[end] && error("axial sample outside pressure-centre span")
    j==nz && return (1-alpha)*matrix[i,j]+alpha*matrix[next,j]
    beta=(query_z-z[j])/(z[j+1]-z[j])
    (1-alpha)*((1-beta)*matrix[i,j]+beta*matrix[i,j+1])+
        alpha*((1-beta)*matrix[next,j]+beta*matrix[next,j+1])
end
function coordinate_mismatch(coarse,fine)
    a=coarse.a;fa=fine.a
    @assert fa.grid.n_theta==2a.grid.n_theta && fa.grid.n_axial==2a.grid.n_axial
    sampled=[sample_pressure(fa.axes.theta_centres,fa.axes.z_centres,fine.pressure,t,z)
        for t in a.axes.theta_centres,z in a.axes.z_centres]
    weights=a.grid.areas
    sqrt(sum((coarse.pressure-vec(sampled)).^2 .* weights)/sum(vec(sampled).^2 .* weights))
end
