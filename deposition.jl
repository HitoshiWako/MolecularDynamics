
using Pkg
using Plots
using Printf

Pkg.activate("MolecularDynamics")
using MolecularDynamics

kB = 1.380649e-2 # 10-21 J/K

σ = 0.340 # 10-9 m [nm]
ϵ = 1.67 # 10-21 J 
m = 66.34 # 10-27 kg

λ = 0.4

Δt = 0.01

n = 12000
samp = 10
interval = 1000

row = 3
col = 15

v = 400

cos_rand() = asin(2*rand()-1)

ps = []
append!(ps,vec([Particle([(x-(col+1)/2)*λ/σ,-(row-2y)*√3λ/2σ],  [0.0,0.0],1.0) for x=1:col,y=1:div(row  ,2)]))
append!(ps,vec([Particle([(x-col/2)*λ/σ,    -(row-2y+1)*√3λ/2σ],[0.0,0.0],1.0) for x=1:col,y=1:div(row+1,2)]))
fixs = [Particle([(x-(col+1)/2)*λ/σ,-row*√3λ/2σ],[0.0,0.0],1.0) for x=1:col]

anim =Animation()
@time for i=0:n
    if i%interval == 0
        θ = cos_rand()-π/2
        push!(ps,Particle([(rand()-0.5)*λ*col/σ,10λ/σ],[v*√(m/ϵ)*1e-3*cos(θ),v*√(m/ϵ)*1e-3*sin(θ)],1.0))
    end   
    if i%samp == 0
        xs = [p.x[1]*σ for p in ps]
        ys = [p.x[2]*σ for p in ps]
        plt = plot(xs,ys,st=:scatter,label="",title=@sprintf("%3.1f ps",i*Δt*√(m*σ^2/ϵ)),
                    xlabel="nm",ylabel="nm",xlims=(-λ*col/2,λ*col/2),ylims=(-λ*(row+1),12λ),aspect_ratio=:equal)
        fxs = [p.x[1]*σ for p in fixs]
        fys = [p.x[2]*σ for p in fixs]
        plot!(fxs,fys,st=:scatter,label="")

        frame(anim,plt)
    end

    accs = calc_acc(lennardjones_force,ps,fixs,[-col*λ/2σ,col*λ/2σ])
    update_x!(ps,accs,Δt)
    fold!(ps,[-col*λ/2σ,col*λ/2σ])
    accs_next = calc_acc(lennardjones_force,ps,fixs,[-col*λ/2σ,col*λ/2σ])
    update_v!(ps,accs,accs_next,Δt)

end
gif(anim,"results/deposition_3x15.gif")