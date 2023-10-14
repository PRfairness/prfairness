using MatrixNetworks

struct Graph_direct_
    n :: Int # |V|
    m :: Int # |E|
    u :: Array{Int, 1}
    v :: Array{Int, 1} # uv is an edge
#    w :: Array{Float64, 1} # weight of each edge
    nbr_in  :: Array{Array{Int, 1}, 1}
    nbr_out :: Array{Array{Int, 1}, 1}
    degree::Vector{Int}
    S::Vector{Int}
    name :: String
end

include("core.jl")

function get_graph_direct(ffname)
    n = 0
    n2 = 0
    Label = Dict{Int32, Int32}()
    Label2 = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = n += 1
    getID2(x :: Int) = haskey(Label2, x) ? Label2[x] : Label2[x] = n2 += 1

    fname = string("./graphs/",ffname,".txt")
    fin = open(fname, "r")

    str = readline(fin)
    u = Int[]
    v = Int[]
    str = strip(str)
    tot = 0
    i=0
    while str != ""
        if occursin("#", str)
            str = strip(readline(fin))
            continue
        end
        str = split(str)
        if length(str)!=2
            str = strip(readline(fin))
            continue
        end
        i = i + 1
        x   = parse(Int, str[1])
        y   = parse(Int, str[2])
        if x!=y
            u1 = getID(x)
            v1 = getID(y)
            Origin[u1] = x
            Origin[v1] = y
            push!(u, u1)
            push!(v, v1)
            tot += 1
        end
        str = strip(readline(fin))
    end
    nbr_in=[ [ ] for i in 1:n ]
    nbr_out=[ [ ] for i in 1:n ]
    for i=1:tot
        u1=u[i];
        v1=v[i];
        push!(nbr_out[u1],v1);
        push!(nbr_in[v1],u1);
    end

    close(fin)

    # our_degree = zeros(n)
    # zeros_out = []
    # for i=1:n
    #     our_degree[i] = length(nbr_out[i])
    #     if our_degree[i] == 0
    #         push!(zeros_out, i)
    #     end
    # end
    # println(n, " ", tot)

    # #remove out degree == 0
    # # n2 = n - length(zeros_out)
    # tot2   = 0
    # u2 = Int[]
    # v2 = Int[]

    # for i=1:tot
    #     u1=u[i];
    #     v1=v[i];
    #     if !(u1 in zeros_out || v1 in  zeros_out)
    #         u1 = getID2(u1)
    #         v1 = getID2(v1)
    #         push!(u2, u1)
    #         push!(v2, v1)
    #         tot2 += 1
    #     end
    # end
    # nbr_in2=[ [ ] for i in 1:n2 ]
    # nbr_out2=[ [ ] for i in 1:n2 ]
    # for i=1:tot2
    #     u1=u2[i];
    #     v1=v2[i];
    #     push!(nbr_out2[u1],v1);
    #     push!(nbr_in2[v1],u1);
    # end

    # our_degree2 = zeros(n2)
    # for i=1:n2
    #     our_degree2[i] = length(nbr_out2[i])
    # end
    #read attribute

    fname = string("./graphs/",ffname,"_c.txt")
    fin = open(fname, "r")
    str = strip(readline(fin))
    S = []
    iv = 0
    while str != ""
        # if ! occursin(" ", str) & (! occursin("\t", str))
        #     str = strip(readline(fin))
        #     continue
        # end
        
        iv = iv + 1
        if iv % 10000 == 0 
            println(str, " ", iv)
        end
        str = split(str)

        if length(str) == 2
            x   = parse(Int, str[1])
            y   = parse(Int, str[2])
            # if  haskey(Label2, x) && y == 0
            #     push!(S, Label2[x])
            # end
            if  haskey(Label, x) && y == 0
                    push!(S, Label[x])
            end
        end
        str = strip(readline(fin))
    end
    close(fin)

    return Graph_direct_(n, tot, u, v, nbr_in, nbr_out, our_degree, S, ffname)
end

function remove_edge!(g, u, v)
    index = findfirst(isequal(v), g.nbr_out[u])
    splice!(g.nbr_out[u], index)

    index = findfirst(isequal(u), g.nbr_in[v])
    splice!(g.nbr_in[v], index)

    index = findfirst((g.u .== u) .& (g.v .== v))
    splice!(g.u, index)
    splice!(g.v, index)
end

function add_edge!(g, u, v)
    push!(g.nbr_out[u], v)
    push!(g.nbr_in[v], u)
    push!(g.u, u)
    push!(g.v, v)
end
