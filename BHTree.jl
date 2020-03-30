module BHTree

using Bodies, LinearAlgebra

export OctTree, populate!, tree_evolve!, maxdepth, mindepth

mutable struct OctTree
    size::Float64
    origin::Vector{Float64}
    body::Union{Body, Nothing}
    center::Center
    children::Vector{OctTree}

    OctTree(s::Float64, o::Vector{Float64}, c::Bool) = new(s, o, nothing, Center(o, 0., c), Vector{OctTree}())
    OctTree(s::Float64, o::Vector{Float64}) = new(s, o, nothing, Center(o, 0., true), Vector{OctTree}())
    OctTree(s::Float64, c::Bool) = new(s, zeros(3), nothing, Center(zeros(3), 0., c), Vector{OctTree}())
end

isbody(t::OctTree) = t.body != nothing
isempty(t::OctTree) = t.center.m == 0.

maxdepth(t::OctTree) = isempty(t) ? 0 : (isbody(t) ? 1 : 1 + maximum([maxdepth(c) for c in t.children]))
mindepth(t::OctTree) = isempty(t) ? 0 : (isbody(t) ? 1 : 1 + minimum([mindepth(c) for c in t.children]))

populate!(t::OctTree, bs::Vector{Body}) = for b in bs populate!(t, b) end

function populate!(t::OctTree, b::Body)
    if isempty(t)
        t.body = b
    elseif isbody(t)
        t.children = children(t.origin, t.size / 2)
        oct1 = octant(t.body.q, t.origin)
        oct2 = octant(b.q, t.origin)
        populate!(t.children[oct1], t.body)
        populate!(t.children[oct2], b)
        t.body = nothing
    else
        oct = octant(b.q, t.origin)
        populate!(t.children[oct], b)
    end
    add!(t.center, b)
end

function children(o::Vector{Float64}, s::Float64)
    Q = hcat(o .+ s/2, o .- s/2)
    vec([OctTree(s, [x, y, z]) for x in Q[1,:], y in Q[2,:], z in Q[3,:]])
end

function forces(t::OctTree, g::OctTree, θ::Float64)
    Fs = []
    if isbody(t)
        push!(Fs, (force(t.body, g, θ), t.body.i))
    else
        cs = filter(!(isempty), t.children)
        Threads.@threads for c in cs Fs = vcat(Fs, forces(c, g, θ)) end
    end
    Fs
end

function force(b::Body, t::OctTree, θ::Float64)
    F = zeros(3)
    if !(isbody(t) && b.i == t.body.i)
        s = t.size
        d = norm(t.center.q - b.q)
        if θ < s/d
            if isbody(t)
                F = gravity(b, t.body)
            else
                cs = filter(!(isempty), t.children)
                Threads.@threads for c in cs F += force(b, c, θ) end
            end
        else
            F = gravity(b, t.center)
        end
    end
    F
end

function tree_evolve!(bs::Vector{Body}, dt::Float64, L::Float64, θ::Float64)
    for b in bs
        b.q .+= b.p / b.m * 0.5dt
    end
    tree = OctTree(L, true)
    populate!(tree, bs)
    Fs = forces(tree, tree, θ)
    for (F, i) in Fs
        bs[i].p .+= F * dt
    end
    for b in bs
        b.q .+= b.p / b.m * 0.5dt
    end
end

function octant(q, origin)
    x = q - origin
    if sign(x[3]) == 1
        if sign(x[2]) == 1
            return sign(x[1]) == 1 ? 1 : 2
        else
            return sign(x[1]) == 1 ? 3 : 4
        end
    else
        if sign(x[2]) == 1
            return sign(x[1]) == 1 ? 5 : 6
        else
            return sign(x[1]) == 1 ? 7 : 8
        end
    end
end


end
