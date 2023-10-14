using Random
using Dates
using DataStructures
include("logw.jl")
include("graph.jl")

seed_value = Dates.millisecond(now())
Random.seed!(seed_value)


function heapsearch(a, k)

    queue = Queue{Tuple{Float64, Int}}() 
    for (i, val) in enumerate(a)
        if length(queue) < k || val > first(queue)[1]
            if length(queue) == k
                dequeue!(queue)
            end
            enqueue!(queue, (val, i))
        end
    end
    return queue
end

function loop_erased_walk(g::Graph_direct_, alpha)
    n = g.n
    intree = falses(n)
    next = fill(-1, n)
    root = fill(-1, n)
    for i in 1:n
        if g.degree[i] == 0
            continue
        end
        u = i
        while !intree[u]
            if rand() < 1/(1 + (1 - alpha)/alpha)
                intree[u] = true
                root[u] = u
                next[u] = -1
            else
                next[u] = g.nbr_out[u][rand(1:g.degree[u])]
                u = next[u]
            end
        end
        r = root[u]
        u = i
        while !intree[u]
            root[u] = r
            intree[u] = true
            u = next[u]
        end
    end
    w = 1
    for i=1:n
        if i != root[i]
            w = w * 1/(alpha/(1 - alpha)*g.degree[i])
        end
    end
    return root, w
end


function fastv(g, v, budget, psi, alpha)
    n = g.n
    m = g.m
    # g = findconnect_direct_strong(g)
    
    S = copy(g.S)
    
    sigma = zeros(g.n)
    tsigma = zeros(g.n)
    eta = zeros(g.n)
    ssss = zeros(Int, g.n)
    for i in 1:length(S)
        ssss[S[i]] = 1
    end
    deg = zeros(n)
    dmax = maximum(deg)
    R = []

    for b=1:budget

        A = adjsp_direct(g);
        rr = 0.0
        dm1 = spdiagm([1 / i for i in g.degree])
        pm = dm1 * A

        for rou in 1:psi
            root. w = loop_erased_walk(g, alpha)
            rr += w
            for j in 1:g.n
                if ssss[u] == 1
                    eta[j] += w
                end
            end
            u = root[v]
            tsigma[u] += w
        end

        for i in 1:length(tsigma)
            tsigma[i] = tsigma[i] / rr
            eta[i] = eta[i] / rr
        end
        kcand = heapsearch(eta, dmax)
        
        delta = 0
        ii = 0; jj = 0; kk = 0;
        idx = 1
        for a in 1:m
            i = g.u[a]
            j = g.v[a]
            for (val, pos) in kcand
                if A[i,pos] == 0
                    c_delta = ((1 - alpha) * pm[i,j] / n  * tsigma[i] * (val - eta[j]))
                    if c_delta > delta
                        idx = a 
                        ii = i; jj = j; kk = pos;
                        delta = c_delta
                    end
                end
            end
        end
        push!(R, (ii, jj, kk))
        remove_edge!(g, ii, jj)
        add_edge!(g, ii, kk)
    end
    return R
end

str = "books"
v = 1
budget = 50
psi = 1000
alpha = 0.1
g = get_graph_direct(str)
fastv(g, v, budget, psi, alpha )