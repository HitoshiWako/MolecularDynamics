using Pkg
using Plots
using Printf

Pkg.activate("MolecularDynamics")
using MolecularDynamics

σ = 0.340
ϵ = 1.67 # 1.67
m = 66.34 # 10-27kg

λ = 0.4

Δt = 0.01

n = 10000
samp = 100

row = 30
col = 30
v = 4.0

ps =[Particle([0.0,5λ/σ],[0.0,-v*√(m/ϵ)],1.0)]
append!(ps,vec([Particle([(x-(col+1)/2)*λ/σ,-(row-2y)*√3λ/2σ],  [0.0,0.0],1.0) for x=1:col,y=1:div(row  ,2)]))
append!(ps,vec([Particle([(x-col/2)*λ/σ,    -(row-2y+1)*√3λ/2σ],[0.0,0.0],1.0) for x=1:col,y=1:div(row+1,2)]))
fixs = [Particle([(x-(col+1)/2)*λ/σ,-row*√3λ/2σ],[0.0,0.0],1.0) for x=1:col]


anim =Animation()
@time for i=0:n
    if i%samp == 0
        xs = [p.x[1]*σ for p in ps]
        ys = [p.x[2]*σ for p in ps]
        plt = plot(xs,ys,st=:scatter,label="",title=@sprintf("%3.1f ps",i*Δt*√(m*σ^2/ϵ)),
                    xlabel="nm",ylabel="nm",xlims=(-λ*col/2,λ*col/2),ylims=(-λ*(row+1),6λ),aspect_ratio=:equal)
        fxs = [p.x[1]*σ for p in fixs]
        fys = [p.x[2]*σ for p in fixs]
        plot!(fxs,fys,st=:scatter,label="")
        frame(anim,plt)
    end
    accs = calc_acc(lennardjones_force,ps,fixs,[-col*λ/2σ,col*λ/2σ])
    update_x!(ps,accs,Δt)
    fold!(ps,[-col*λ/2σ,col*λ/2σ])
    accs_next =calc_acc(lennardjones_force,ps,fixs,[-col*λ/2σ,col*λ/2σ])
    update_v!(ps,accs,accs_next,Δt)
end
gif(anim,"results/collision_fixed.gif")