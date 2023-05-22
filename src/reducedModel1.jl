
mutable struct ModelR1 <: AbstractMCESModel
    model::Model
    x::Array{VariableRef,2}
    zA::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Tuple{Int64, Int64, Int64}}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64, Int64}, Int64}}}}
    zB::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Tuple{Int64, Int64, Int64}}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64, Int64}, Int64}}}}
    # zB::Dict{Array{Int64,1}, VariableRef}
    gA::Graph{Int64}
    gB::Graph{Int64}
end

@kwdef struct ConfigR1 <: Config
    name::String
    optimizer::String = "Gurobi"
    max_runtime::Int = 1000
    relaxation::Bool=false
    formulation::String = "R1"
    adjacency2::Bool = false
    adjacency3::Bool = false
end

function add_adjacency(m::ModelR1)
    for u2 in vertices(m.gA), v1 in vertices(m.gB)
        @constraint(m.model, sum(m.zA[(u1,u2,v1)] for u1 in inneighbors(m.gA, u2)) <= sum(m.zB[(v2, v1, u2)] for v2 in outneighbors(m.gB, v1)) )
    end

    for u1 in vertices(m.gA), u2 in outneighbors(m.gA, u1), v1 in vertices(m.gB)
        @constraint(m.model, m.zA[(u1, u2, v1)] <= m.x[u1, v1])
    end
    for v1 in vertices(m.gB), v2 in outneighbors(m.gB, v1), u1 in vertices(m.gA)
        @constraint(m.model, m.zB[(v1, v2, u1)] <= m.x[u1, v1])
    end
end

function add_adjacency2(m::ModelR1)
    for u1 in vertices(m.gA), v1 in vertices(m.gB)
        @constraint(m.model, sum(m.zA[(u1,u2,v1)] for u2 in neighbors(m.gA, u1)) == sum(m.zB[(v1, v2, u1)] for v2 in neighbors(m.gB, v1)))
    end
end

function add_adjacency3(m::ModelR1)
    for u1 in vertices(m.gA), u2 in neighbors(m.gA, u1)
        u1 <= u2 && continue
        @constraint(m.model, sum(m.zA[(u1, u2, v1)] for v1 in vertices(m.gB)) == sum(m.zA[(u2, u1, v1)] for v1 in vertices(m.gB)))
    end
    for v1 in vertices(m.gB), v2 in neighbors(m.gB, v1)
        v1 <= v2 && continue
        @constraint(m.model, sum(m.zB[(v1, v2, u1)] for u1 in vertices(m.gA)) == sum(m.zB[(v2, v1, u1)] for u1 in vertices(m.gA)))
    end
end

function run_config(cf::ConfigR1, gA::AbstractGraph, gB::AbstractGraph)
    model = get_model(cf)

    clean_graph!(gA)
    clean_graph!(gB)
    # gA must be the smallest graph
    if nv(gB) < nv(gA)
        gA, gB = gB, gA
    end

    zA_indices = [(src(e), dst(e), v) for e in all_edges(gA), v in vertices(gB) if src(e) != dst(e)]
    zB_indices = [(src(e), dst(e), u) for e in all_edges(gB), u in vertices(gA) if src(e) != dst(e)]
    if cf.relaxation
        @variable(model, 0<=x[i=vertices(gA),j=vertices(gB)]<=1)
        @variable(model, 0<=zA[zA_indices])
        @variable(model, 0<=zB[zB_indices])
    else
        @variable(model, x[i=vertices(gA),j=vertices(gB)], Bin)
        @variable(model, zA[zA_indices], Bin)
        @variable(model, zB[zB_indices], Bin)
    end

    m = ModelR1(model,x, zA, zB, gA, gB)
    loop_cons = handle_loops!(gA, gB, x)
    obj = sum(m.zA[(e.src, e.dst, v)] + m.zA[(e.dst, e.src, v)] for e in edges(m.gA), v in vertices(m.gB)) + 2*loop_cons
    @objective(m.model, Max, obj)

    add_injection(m)
    add_adjacency(m)
    cf.adjacency2 && add_adjacency2(m)
    cf.adjacency3 && add_adjacency3(m)

    optimize!(m.model)
    display_solution(m, 2)

    return m, 2
end
