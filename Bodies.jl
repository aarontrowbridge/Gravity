module Bodies

using LinearAlgebra

export Body, Center, gravity, add!, erase!, collisions!, trim!, coalesce!, coalesce_star!, prep!, brute_evolve!

# body struct
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

# center of mass struct
mutable struct Center
    q::Vector{Float64}
    m::Float64
    c::Bool
end

# radius function R : mass -> radius
R(m) = (3m / (4π))^(1/3)

# gravitational constant
const G = 4π^2

# gravitational damping constant
const A = 0.0001

# short range repulsion constant
const β = 1.0

# magnitude of gravity
mag(m1, m2, r) = G * m1 * m2 / (r^2 + A^2) - β * r^(-6)

# force of gravity on body a from body b
function gravity(a::Body, b::Body)
    r = norm(b.q - a.q)
    u = (b.q - a.q) / r
    F = mag(a.m, b.m, r) * u
    if xor(a.c, b.c)
        return a.c ? -F : F
    elseif a.c & b.c
        return F
    else
        return -F
    end
end


# force of gravity on body b from center of mass c
function gravity(b::Body, c::Center)
    r = norm(c.q - b.q)
    u = (c.q - b.q) / r
    F = mag(b.m, c.m, r) * u
    if xor(b.c, c.c)
        return b.c ? -F : F
    elseif b.c & c.c
        return F
    else
        return -F
    end
end

# add body b & a to get new body
Base.:+(a::Body, b::Body) = Body((a.m * a.q + b.m * b.q) / (a.m + b.m), a.p + b.p, a.m + b.m)

# add body b to body a
function add!(a::Body, b::Body)
    m = a.m + b.m
    a.q = (a.m * a.q + b.m * b.q) / m
    a.p += b.p
    a.m = m
    a.r = R(m)
    erase!(b)
end

# add body b to center c
add!(c::Center, b::Body) = begin M = c.m + b.m; c.q = (c.m * c.q + b.m * b.q) / M; c.m = M end

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
            if R < a.r + b.r && !(pair in pairs)
                push!(pairs, pair)
                add!(a, b)
            end
        end
    end
    filter!(b -> b.i != 0, bs)
    for k = 1:length(bs) bs[k].i = k end
end

function coalesce_star!(bs::Vector{Body})
    for j = 1:length(bs)-1
        if norm(bs[j].q - bs[end].q) < bs[j].r + bs[end].r
            add!(bs[end], bs[j])
        end
    end
    filter!(b -> b.i != 0, bs)
    for k = 1:length(bs) bs[k].i = k end
end

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
