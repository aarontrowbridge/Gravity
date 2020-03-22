push!(LOAD_PATH, pwd())

using Bodies, BHTree, LinearAlgebra, Distributions

function spawn(c)
    bs = Vector{Body}()
    for i = 1:N
        x = rand(Normal(0, L/8)) # rand(-L/2:.001:L/2)
        y = rand(Normal(0, L/8)) # rand(-L/2:.001:L/2)
        z = 0
        m = m₀
        q = [x, y, z]
        p = m * cross(q, [0, 0, ω])
        b = Body(q, p, m, c, i)
        push!(bs, b)
    end
    bs
end

const N = 500
const L = 200.
const m₀ = 0.1
const ω = 0.1
const dt = 0.001
const θ = 1.0
const star = false
const maxitr = 1e5

function main()
    bodies::Vector{Body} = spawn(true)
    prep!(bodies, L)

    if star
        sun = Body(zeros(3), zeros(3), 100., length(bodies) + 1)
        push!(bodies, sun)
    end

    frame = 0
    while frame <= maxitr

        # brute_evolve!(bodies, dt)
        tree_evolve!(bodies, dt, L, θ)

        trim!(bodies, L)

        coalesce!(bodies)

        if frame % 50 == 0
            anim.(bodies)
            println("F")
        end
        frame += 1

    end
end

main()

