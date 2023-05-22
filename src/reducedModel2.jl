
mutable struct ModelR2 <: AbstractMCESModel
    model::Model
    x::Array{VariableRef,2}
    z::Array{VariableRef,2}
    gA::Graph{Int64}
    gB::Graph{Int64}
end

@kwdef struct ConfigR2 <: Config
    name::String
    optimizer::String = "Gurobi"
    max_runtime::Int = 1000
    relaxation::Bool = false
    add_degree_cuts::Bool = false
    formulation::String = "R2"
end

function add_adjacency(m::ModelR2)
    for u2 in vertices(m.gA), v1 in vertices(m.gB)
        @constraint(m.model, m.z[u2, v1] <= sum(m.x[u1, v1] for u1 in inneighbors(m.gA, u2)) )
        if is_directed(m.gA)
            @constraint(m.model, m.z[u2, v1] <= sum(m.x[u2, v2] for v2 in outneighbors(m.gB, v1)) )
        end
    end
    for u1 in vertices(m.gA), v2 in vertices(m.gB)
        @constraint(m.model, m.z[u1, v2] <= sum(m.x[u1, v1] for v1 in inneighbors(m.gB, v2)) )
        if is_directed(m.gA)
            @constraint(m.model, m.z[u1, v2] <= sum(m.x[u2, v2] for u2 in outneighbors(m.gA, u1)) )
        end
    end
end

function add_degree_cuts(m::ModelR2)
    for u in vertices(m.gA)
        @constraint(m.model, sum(m.z[u,v] for v in vertices(m.gB)) <= sum(m.x[u,v]*min(degree(m.gA, u), degree(m.gB, v)) for v in vertices(m.gB)) )
    end
    for v in vertices(m.gB)
        @constraint(m.model, sum(m.z[u,v] for u in vertices(m.gA)) <= sum(m.x[u,v]*min(degree(m.gA, u), degree(m.gB, v)) for u in vertices(m.gA)) )
    end
end

function run_config(cf::ConfigR2, gA::AbstractGraph, gB::AbstractGraph)
    model = get_model(cf)

    clean_graph!(gA)
    clean_graph!(gB)
    # gA must be the smallest graph
    if nv(gB) < nv(gA)
        gA, gB = gB, gA
    end

    if cf.relaxation
        @variable(model, 0<=x[i=vertices(gA),j=vertices(gB)]<=1)
        @variable(model, 0<=z[i=vertices(gA),j=vertices(gB)])
    else
        @variable(model, x[i=vertices(gA),j=vertices(gB)], Bin)
        @variable(model, z[i=vertices(gA),j=vertices(gB)], Bin)
    end

    m = ModelR2(model, x, z, gA, gB)
    loop_cons = handle_loops!(gA, gB, x)
    obj = sum(m.z[u, v] for u in vertices(m.gA), v in vertices(m.gB)) + 2*loop_cons
    @objective(m.model, Max, obj)

    add_injection(m)
    add_adjacency(m)
    cf.add_degree_cuts && add_degree_cuts(m)

    optimize!(m.model)
    display_solution(m, 2)

    return m, 2
end
