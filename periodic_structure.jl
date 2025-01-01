using Pkg
using Plots
using Printf
using LinearAlgebra

Pkg.activate("MolecularDynamics")
using MolecularDynamics

function calc_acc(func,ps,fixs,bound,cutoff=8)
    alls =vcat(ps,fixs)
    accs = []
    for p1 in ps
        acc = zeros(size(p1.x))
        for  p2 in alls
            if p1 != p2 
                if norm(p1.x-p2.x) < cutoff
                    acc += func(norm(p1.x-p2.x))/p1.m*(p1.x-p2.x)/norm(p1.x-p2.x)
                elseif norm(p1.x-p2.x-[bound,0.0]) < cutoff
                    acc += func(norm(p1.x-p2.x-[bound,0.0]))/p1.m*(p1.x-p2.x-[bound,0.0])/norm(p1.x-p2.x-[bound,0.0])
                elseif norm(p1.x-p2.x+[bound,0.0]) < cutoff
                    acc += func(norm(p1.x-p2.x+[bound,0.0]))/p1.m*(p1.x-p2.x+[bound,0.0])/norm(p1.x-p2.x+[bound,0.0])
                end
            end
        end
        push!(accs,acc)
    end
    return accs
end

function fold!(ps,bound)
    for p in ps
        if p.x[1] > bound
            p.x[1] -= bound
        elseif p.x[1] < 0
            p.x[1] += bound
        end
    end
    return ps
end

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
    accs = calc_acc(lennardjones_force,ps,fixs,col*λ/σ,3)
    update_x!(ps,accs,Δt)
    fold!(ps,col*λ/σ)
    accs_next =calc_acc(lennardjones_force,ps,fixs,col*λ/σ,3)
    update_v!(ps,accs,accs_next,Δt)
end
gif(anim,"results/periodic.gif")