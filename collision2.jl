using Pkg
using Plots
using Printf
using LinearAlgebra

Pkg.activate("MolecularDynamics")
using MolecularDynamics

function calc_acc(func,ps,fixs,bound=(0.0,0.0),cutoff=3)
    alls =vcat(ps,fixs)
    accs = []
    for p1 in ps
        acc = zeros(size(p1.x))
        for  p2 in alls
            if p1 != p2 
                if norm(p1.x-p2.x) < cutoff
                    acc += func(norm(p1.x-p2.x))/p1.m*(p1.x-p2.x)/norm(p1.x-p2.x)
                elseif norm(p1.x-p2.x-[bound[2]-bound[1],0.0]) < cutoff
                    acc += func(norm(p1.x-p2.x-[bound[2]-bound[1],0.0]))/p1.m*(p1.x-p2.x-[bound[2]-bound[1],0.0])/norm(p1.x-p2.x-[bound[2]-bound[1],0.0])
                elseif norm(p1.x-p2.x+[bound[2]-bound[1],0.0]) < cutoff
                    acc += func(norm(p1.x-p2.x+[bound[2]-bound[1],0.0]))/p1.m*(p1.x-p2.x+[bound[2]-bound[1],0.0])/norm(p1.x-p2.x+[bound[2]-bound[1],0.0])
                end
            end
        end
        push!(accs,acc)
    end
    return accs
end

function fold!(ps,bound)
    for p in ps
        if p.x[1] > bound[2]
            p.x[1] -= bound[2]-bound[1]
        elseif p.x[1] < bound[1]
            p.x[1] += bound[2]-bound[1]
        end
    end
    return ps
end

σ = .340
ϵ = 1.67
m = 66.34 # 10-27kg
v = 3.0

λ = 0.4

Δt = 0.01

n = 1000
samp = 10

row = 40
col = 40

ps =[Particle([0.0,λ/σ],[0.0,-v*√(m/ϵ)],1.0)]
append!(ps,vec([Particle([(x-(col+1)/2)*λ/σ,-(row-2y)*√3λ/2σ],  [0.0,0.0],1.0) for x=1:col,y=1:div(row  ,2)]))
append!(ps,vec([Particle([(x-col/2)*λ/σ,    -(row-2y+1)*√3λ/2σ],[0.0,0.0],1.0) for x=1:col,y=1:div(row+1,2)]))
fixs = [Particle([(x-(col+1)/2)*λ/σ,-row*√3λ/2σ],[0.0,0.0],1.0) for x=1:col]


anim =Animation()
@time for i=0:n
    if i%samp == 0
        xs = [p.x[1]*σ for p in ps]
        ys = [p.x[2]*σ for p in ps]
        plt = plot(xs,ys,st=:scatter,label="",title=@sprintf("%3.1f ps",i*Δt*√(m*σ^2/ϵ)),
                    xlabel="nm",ylabel="nm",xlims=(-λ*col/2,λ*col/2),ylims=(-λ*(row+1),5λ),aspect_ratio=:equal)
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
gif(anim,"results/collision2.gif")