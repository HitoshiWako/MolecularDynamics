using Pkg
using Plots
using Printf

Pkg.activate("MolecularDynamics")
using MolecularDynamics

σ = 0.340
ϵ = 1.67
m = 66.34 # 10-27kg

λ = 0.4

Δt = 0.01

n = 10000
samp = 10

row = 10
col = 10
ps = vec([Particle([(x-1)*λ/σ,(y-1)*λ/σ],[0.0,0.0],1.0) for x=1:col,y=1:row])

anim =Animation()
@time for i=0:n
    if i%samp == 0
        xs = [p.x[1]*σ for p in ps]
        ys = [p.x[2]*σ for p in ps]
        plt = plot(xs,ys,st=:scatter,label="",title=@sprintf("%3.1f ps",i*Δt*√(m*σ^2/ϵ)),
                    xlabel="nm",ylabel="nm",xlims=(-λ,λ*col),ylims=(-λ,λ*row),aspect_ratio=:equal)
        frame(anim,plt)
    end
    accs = calc_acc(lennardjones_force,ps)
    update_x!(ps,accs,Δt)
    accs_next =calc_acc(lennardjones_force,ps)
    update_v!(ps,accs,accs_next,Δt)
end
gif(anim,"results/lattice.gif")