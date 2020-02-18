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

function brute_evolve!(bs::Vector{Body})
    for b in bs
        b.q .+= b.p / b.m * 0.5dt
    end
    collisions!(bs)
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
    collisions!(bs)
end

function tree_evolve!(bs::Vector{Body})
    for b in bs
        b.q .+= b.p / b.m * 0.5dt
    end
    collisions!(bs)
    tree = Tree(L, true)
    populate!(tree, bs)
    Fs = forces(tree, tree, θ)
    for (F, i) in Fs
        bs[i].p .+= F * dt
    end
    for b in bs
        b.q .+= b.p / b.m * 0.5dt
    end
    collisions!(bs)
end

const N = 100
const L = 200.
const m₀ = 0.1
const ω = 0.1
const dt = 0.001
const θ = 0.50
const star = true
const maxitr = 1e6

function main()
    bodies::Vector{Body} = spawn(true)
    prep!(bodies, L)

    if star
        sun = Body(zeros(3), zeros(3), 100., length(bodies) + 1)
        push!(bodies, sun)
    end

    # io = open("$(pwd())/data/n$(N)_$(star ? star : nostar)_theta$(θ).dat", "w")

    frame = 0
    while frame <= maxitr
        # brute_evolve!(bodies)
        tree_evolve!(bodies)

        trim!(bodies, L)

        if frame % 100 == 0
            anim(bodies)
            println("F")
            # write(io, "F\n")
        end

        frame += 1
    end
end

main()

