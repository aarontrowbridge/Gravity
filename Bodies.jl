module Bodies

using LinearAlgebra, Printf

export Body, anim, gravity, octant, add!, erase!, collisions!, trim!, coalesce!, prep!, brute_evolve!

mutable struct Body
    q::Vector{Float64}
    p::Vector{Float64}
    m::Float64
    r::Float64
    c::Bool
    i::Int64

    Body(q::Vector{Float64},
         p::Vector{Float64},
         m::Float64,
         c::Bool,
         i::Int64) = new(q, p, m, R(m), c, i)

    Body(q::Vector{Float64},
         p::Vector{Float64},
         m::Float64,
         c::Bool) = new(q, p, m, R(m), c, 0)

    Body(q::Vector{Float64},
         p::Vector{Float64},
         m::Float64,
         i::Int64) = new(q, p, m, R(m), true, i)

    Body(q::Vector{Float64},
         p::Vector{Float64},
         m::Float64) = new(q, p, m, R(m), true, 0)

    Body(q::Vector{Float64},
         m::Float64,
         c::Bool) = new(q, zeros(3), m, R(m), c, 0)
end

R(m) = (3m / (4π))^(1/3)

const G = 4π^2

function gravity(a::Body, b::Body)
    r = norm(b.q - a.q)
    u = (b.q - a.q) / r
    mag = G * a.m * b.m / r^2
    F = mag * u
    if xor(a.c, b.c)
        return a.c ? -F : F
    elseif a.c & b.c
        return F
    else
        return -F
    end
end

function gravity(b::Body, x::Vector{Float64}, m::Float64, c::Bool)
    r = norm(x - b.q)
    u = (x - b.q) / r
    mag = G * b.m * m / r^2
    F = mag * u
    if xor(b.c, c)
        return b.c ? -F : F
    elseif b.c & c
        return F
    else
        return -F
    end
end

Base.:+(a::Body, b::Body) = Body((a.m * a.q + b.m * b.q) / (a.m + b.m), a.p + b.p, a.m + b.m)

function add!(a::Body, b::Body)
    m = a.m + b.m
    a.q = (a.m * a.q + b.m * b.q) / m
    a.p += b.p
    a.m = m
    a.r = R(m)
    erase!(b)
end

function erase!(b::Body)
    b.q = zeros(3)
    b.p = zeros(3)
    b.m = 0.
    b.r = 0.
    b.i = 0
end

function collisions!(bs::Vector{Body})
    pairs = []
    for a in bs
        for b in bs
            if a.i == b.i continue end
            pair = Set([a.i, b.i])
            if norm(a.q - b.q) < R(a.m) + R(b.m) && !(pair in pairs)
                push!(pairs, pair)
                μ = 2a.m * b.m / (a.m + b.m)
                x1 = a.q
                x2 = b.q
                v1 = a.p / a.m
                v2 = b.p / b.m
                dsqr = norm(x1 - x2)^2
                a.p -= μ * dot(v1 - v2, x1 - x2) / dsqr * (x1 - x2)
                a.p -= μ * dot(v2 - v1, x2 - x1) / dsqr * (x2 - x1)
            end
        end
    end
end

prep!(bs::Vector{Body}, L::Float64) = begin trim!(bs, L); coalesce!(bs) end

function trim!(bs::Vector{Body}, L)
    filter!(b -> inside(b, L), bs)
    for k = 1:length(bs) bs[k].i = k end
end

inside(body::Body, L) = for q in abs.(body.q) return q < L/2 ? true : false; break end

function coalesce!(bs::Vector{Body})
    pairs = []
    for a in bs
        if a.i == 0 continue end
        for b in bs
            if a.i == b.i || b.i == 0 continue end
            pair = Set([a.i, b.i])
            R = norm(a.q - b.q)
            if R < a.r + b.r && !(pair in pairs) push!(pairs, pair); add!(a, b) end
        end
    end
    filter!(b -> b.i != 0, bs)
    for k = 1:length(bs) bs[k].i = k end
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

anim(b::Body) = b.i == 0 ? nothing : @printf "c3 %f %f %f %f\n" b.q[1] b.q[2] b.q[3] b.r

function brute_evolve!(bs::Vector{Body}, dt::Float64)
    for b in bs
        b.q .+= b.p / b.m * 0.5dt
    end
    for b1 in bs
        F = zeros(3)
        for b2 in bs
            if b2.i != b1.i
                F .+= gravity(b1, b2)
            end
        end
        b1.p .+= F * dt
    end
    for b in bs
        b.q .+= b.p / b.m * 0.5dt
    end
end

end
