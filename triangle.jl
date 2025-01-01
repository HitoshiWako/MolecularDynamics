using Pkg
using Plots
using Printf

Pkg.activate("MolecularDynamics")
using MolecularDynamics

σ = .340
ϵ = 1.67
m = 66.34 # 10-27kg

λ = 0.6

Δt = 0.01

n = 1000
samp = 10

ps = [Particle([-λ/(2σ),-√3λ/(4σ)],[0.0,0.0],1.0),
      Particle([ λ/(2σ),-√3λ/(4σ)],[0.0,0.0],1.0),
      Particle([    0.0, √3λ/(4σ)],[0.0,0.0],1.0)]


anim =Animation()
for i=0:n
    if i%samp == 0
        xs = [p.x[1]*σ for p in ps]
        ys = [p.x[2]*σ for p in ps]
        plt = plot(xs,ys,st=:scatter,label="",title=@sprintf("%3.1f ps",i*Δt*√(m*σ^2/ϵ)),
                    xlabel="nm",ylabel="nm",xlims=(-λ,λ),ylims=(-λ,λ),aspect_ratio=:equal)
        frame(anim,plt)
    end
    accs = calc_acc(lennardjones_force,ps)
    update_x!(ps,accs,Δt)
    accs_next =calc_acc(lennardjones_force,ps)
    update_v!(ps,accs,accs_next,Δt)
end
gif(anim,"results/triangle.gif")