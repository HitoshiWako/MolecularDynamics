module MolecularDynamics
export Particle,lennardjones_potential,lennardjones_force,calc_acc,update_x!,update_v!,fold!

using LinearAlgebra

mutable struct Particle
    x::Vector
    v::Vector
    m
end

lennardjones_potential(r,p=12,q=6)=4*((1/r)^p-(1/r)^q)
lennardjones_force(r,p=12,q=6)=4*(p*(1/r)^(p+1)-q*(1/r)^(q+1))

function calc_acc(func,ps,fixs=[],bound=(0.0,0.0),cutoff=3)
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

function update_x!(ps,accs,Δt)
    for (p,acc) in zip(ps,accs)
        p.x += p.v*Δt+acc/2*Δt^2
    end
    return ps
end

function update_v!(ps,accs,accs_next,Δt)
    for (p,acc,acc_next) in zip(ps,accs,accs_next)
        p.v += (acc+acc_next)/2*Δt
    end
    return ps
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

end # module MolecularDynamics
