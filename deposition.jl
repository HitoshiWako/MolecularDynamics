using Pkg
using Plots
using Printf
using LinearAlgebra

Pkg.activate("MolecularDynamics")
using MolecularDynamics

kB = 1.380649e-2

σ = 0.340
ϵ = 1.67 # 1.67
m = 66.34 # 10-27kg

λ = 0.4

Δt = 0.01

n = 100000
samp = 100
interval = 1000

row = 0
col = 20
temp_vap = 2
temp_sub = 0.1

cos_rand() = acos(1-2rand())

ps = []
append!(ps,vec([Particle([(x-(col+1)/2)*λ/σ,-(row-2y)*√3λ/2σ],  [0.0,0.0],1.0) for x=1:col,y=1:div(row  ,2)]))
append!(ps,vec([Particle([(x-col/2)*λ/σ,    -(row-2y+1)*√3λ/2σ],[0.0,0.0],1.0) for x=1:col,y=1:div(row+1,2)]))
fixs = [Particle([(x-(col+1)/2)*λ/σ,-row*√3λ/2σ],[0.0,0.0],1.0) for x=1:col]

inital_energy = calc_energy(lennardjones_potential, vcat(ps,fixs),[-col*λ/2σ,col*λ/2σ])
anim =Animation()
@time for i=0:n
    if i%interval == 0
        θ = cos_rand()
        v = √(3kB*temp_vap/ϵ)
        push!(ps,Particle([(rand()-0.5)*λ*col,10λ/σ],[v*√(m/ϵ)*cos(θ),-v*√(m/ϵ)*sin(θ)],1.0))
    end   
    if i%samp == 0
        xs = [p.x[1]*σ for p in ps]
        ys = [p.x[2]*σ for p in ps]
        plt = plot(xs,ys,st=:scatter,label="",title=@sprintf("%3.1f ps",i*Δt*√(m*σ^2/ϵ)),
                    xlabel="nm",ylabel="nm",xlims=(-λ*col/2,λ*col/2),ylims=(-λ*(row+1),15λ),aspect_ratio=:equal)
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
    isothermal!(ps,fixs,kB/ϵ*temp_sub)

end
final_energy = calc_energy(lennardjones_potential, vcat(ps,fixs),[-col*λ/2σ,col*λ/2σ])
@printf("Initial:%f,Final:%f,ΔE:%f\n",inital_energy,final_energy,final_energy-inital_energy)
gif(anim,"results/deposition.gif")