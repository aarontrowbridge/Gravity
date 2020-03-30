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

const N = 10
const L = 200.
const m₀ = 1.
const ω = 0.0
const dt = 0.001
const θ = 1.0
const star = false
const maxitr = 1e3

function main()
    println("prepping bodies..")

    bodies::Vector{Body} = spawn(true)
    prep!(bodies, L)

    if star
        sun = Body(zeros(3), zeros(3), 100., length(bodies) + 1)
        push!(bodies, sun)
    end

    println("beginning loops..")

    itr = 0
    while itr <= maxitr

        # brute_evolve!(bodies, dt)
        tree_evolve!(bodies, dt, L, θ)

        trim!(bodies, L)

        # coalesce!(bodies)

        # io = open("data/N$(N)_m$(m₀)_omega$(ω)$(star ? "_star" : "").dat", "w")
        io = stdout

        if itr % 25 == 0
            if itr % 500 == 0
                println("itr = $(itr)")
            end
            anim(bodies, io)
            # println(io, "F")
            # write(io, "F\n")
        end

        itr += 1

    end
end

main()

