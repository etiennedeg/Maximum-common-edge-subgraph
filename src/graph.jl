"""
Some graph utilities
"""

using Graphs: AbstractSimpleGraph, SimpleGraphs.fadj

"""
Bidirectional Edge iterator
For a directed graph, this is a normal edge iterator
For an undirected graph, return all edges twice (once for each orientation)
"""
struct BidirectionalEdgeIter{G} <: AbstractEdgeIter
    g::G
end

Base.eltype(::Type{BidirectionalEdgeIter{SimpleGraph{T}}}) where {T} = Edge{T}
Base.eltype(::Type{BidirectionalEdgeIter{SimpleDiGraph{T}}}) where {T} = Edge{T}

@inline function Base.iterate(eit::BidirectionalEdgeIter{G}, state=(one(eltype(eit.g)), 1) ) where {G <: AbstractSimpleGraph}
    g = eit.g
    fadjlist = fadj(g)
    T = eltype(g)
    n = T(nv(g))
    u, i = state

    n == 0 && return nothing

    @inbounds while true
        list_u = fadjlist[u]
        if i > length(list_u)
            u == n && return nothing

            u += one(u)
            list_u = fadjlist[u]
            i = 1
            continue
        end
        e = Edge(u, list_u[i])
        state = (u, i + 1)
        return e, state
    end

    return nothing
end

Base.length(eit::BidirectionalEdgeIter) = ( is_directed(eit.g) ? ne(eit.g) : 2*ne(eit.g) )

all_edges(g::AbstractSimpleGraph) = BidirectionalEdgeIter(g)

"""
remove unecessary edges (when the formulation allow different graph sizes)
"""
function clean_graph!(graph)
    for u in reverse(vertices(graph))
       if degree(graph, u) == 0
           rem_vertex!(graph, u)
       end
   end
end

"""
remove self-loops (self-loops are reflected in the objective function)
"""
function free_loops!(graph)
    for u in vertices(graph)
        if has_edge(graph, u, u)
            rem_edge!(graph, u, u)
        end
    end
end
