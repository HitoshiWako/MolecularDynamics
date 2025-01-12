using Pkg
using Plots
using Printf
using LinearAlgebra

Pkg.activate("MolecularDynamics")
using MolecularDynamics

σ = .340
ϵ = 1.67
m = 66.34 # 10-27kg

λ = 0.4

Δt = 0.01

n = 10000
samp = 10

row = 10
col = 10

ps = vec([Particle([(x-0.5)*λ/σ,-(y-1)*λ/σ],[0.0,0.0],1.0) for x=1:col,y=1:row])
fixs = [Particle([(x-0.5)*λ/σ,-row*λ/σ],[0.0,0.0],1.0) for x=1:col]


anim =Animation()
@time for i=0:n
    if i%samp == 0
        xs = [p.x[1]*σ for p in ps]
        ys = [p.x[2]*σ for p in ps]
        plt = plot(xs,ys,st=:scatter,label="",title=@sprintf("%3.1f ps",i*Δt*√(m*σ^2/ϵ)),
                    xlabel="nm",ylabel="nm",xlims=(0,λ*col),ylims=(-λ*(row+1),2λ),aspect_ratio=:equal)
        fxs = [p.x[1]*σ for p in fixs]
        fys = [p.x[2]*σ for p in fixs]
        plot!(fxs,fys,st=:scatter,label="")
        frame(anim,plt)
    end
    accs = calc_acc(lennardjones_force,ps,fixs,(0.0,col*λ/σ))
    update_x!(ps,accs,Δt)
    fold!(ps,(0.0,col*λ/σ))
    accs_next =calc_acc(lennardjones_force,ps,fixs,(0.0,col*λ/σ))
    update_v!(ps,accs,accs_next,Δt)
end
gif(anim,"results/periodic.gif")