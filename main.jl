push!(LOAD_PATH, pwd())

using Bodies, BHTree, LinearAlgebra, Distributions

function spawn(c)
    bs = Vector{Body}()
    for i = 1:N
        x = rand(Normal(0, L/8)) # rand(-L/2:0.001:L/2)
        y = rand(Normal(0, L/8)) # rand(-L/2:0.001:L/2)
        z = rand(-L/50:0.001:L/50)
        m = m₀
        q = [x, y, z]
        p = m * cross(q, [0, 0, ω])
        b = Body(q, p, m, c, i)
        push!(bs, b)
    end
    bs
end

function anim(io, bs::Vector{Body})
    for b in bs
        if b.i != 0
            x, y, z = b.q
            r = b.r
            println(io, "c3 $x $y $z $r")
            println(io, "T -0.9 0.9")
            println(io, "n = $(length(bs))")
        end
    end
    println(io, "F")
end

const N = 200
const L = 200.
const m₀ = 0.1
const ω = 0.2
const dt = 0.001
const θ = 2.0
const star = true
const maxitr = 5e5

function main()
    bodies::Vector{Body} = spawn(true)
    prep!(bodies, L)

    if star
        sun = Body(zeros(3), zeros(3), 100., length(bodies) + 1)
        push!(bodies, sun)
    end

    # io = open("data/n$(N)_m$(m₀)_omega$(ω)_theta$(θ)$(star ? "_star" : "").dat", "w")
    io = stdout

    frame = 0
    while frame <= maxitr

        # brute_evolve!(bodies, dt)
        tree_evolve!(bodies, dt, L, θ)

        trim!(bodies, L)

        # coalesce!(bodies)
        coalesce_star!(bodies)

        if frame % 50 == 0
            if frame % 500 == 0 println(frame) end
            anim(io, bodies)
        end

        frame += 1

    end
end

main()

