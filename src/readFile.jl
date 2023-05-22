using Graphs

function load_graph_ARG(path::String, type::Type{T}) where {T<:AbstractGraph}
    open(path) do file
        L = read!(file, Vector{Int16}(undef,stat(file).size ÷ sizeof(UInt16)))
        L = convert.(Int, L)
        n = L[1]
        c = n+2
        succ = Array{Array}(undef,0)
        for i in 1:n
            push!(succ, L[c+1:2:c+2*L[c]].+1)
            c += 2*L[c]+1
        end

        graph = type(n)
        for i in 1:n, j in succ[i]
            add_edge!(graph, i, j)
        end
        return graph
    end
end

function load_graph_ML(path::String, graphType::Type{T}) where {T<:AbstractGraph}
    graph = 0
    # println("salut ", path)
    open(path) do file
        n = parse(Int, readline(file))
        graphA = graphType(n)
        graphB = graphType(n)
        lines = readlines(file)
        filter!(x->!isempty(x), lines)
        succ = map(
            x->map(y->parse(Int,y)+1,x), #convert lines into numerical values
            map(x->split(x), lines) # split lines
        )
        for i in 1:n
            for j in succ[i][2:end]
                add_edge!(graphA, i, j)
            end
        end
        for i in n+1:2*n
            for j in succ[i][2:end]
                add_edge!(graphB, i-n, j)
            end
        end

        return (graphA, graphB)
    end
end
