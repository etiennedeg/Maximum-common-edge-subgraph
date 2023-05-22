
mutable struct ModelD <: AbstractMCESModel
    model::Model
    x::Array{VariableRef,2}
    d::Array{VariableRef,2}
    gA::Graph
    gB::Graph
end

@kwdef struct ConfigD <: Config
    name::String
    optimizer::String = "Gurobi"
    max_runtime::Int = 1000
    relaxation::Bool = false
    formulation::String = "D"
end

function add_adjacency(m::ModelD)
    for u1 in vertices(m.gA), v1 in vertices(m.gB)
        @constraint(m.model, m.d[u1,v1] <= min(degree(m.gA, u1), degree(m.gB, v1))*m.x[u1,v1])
        @constraint(m.model, m.d[u1,v1] <= sum(m.x[u2,v2] for u2 in neighbors(m.gA, u1), v2 in neighbors(m.gB, v1)))
    end
end

function run_config(cf::ConfigD, gA::AbstractGraph, gB::AbstractGraph)
    model = get_model(cf)

    clean_graph!(gA)
    clean_graph!(gB)
    # gA must be the smallest graph
    if nv(gB) < nv(gA)
        gA, gB = gB, gA
    end

    if cf.relaxation
        @variable(model, x[i=vertices(gA),j=vertices(gB)])
        @variable(model, 0<=d[i=vertices(gA),j=vertices(gB)])
    else
        @variable(model, x[i=vertices(gA),j=vertices(gB)], Bin)
        @variable(model, 0<=d[i=vertices(gA),j=vertices(gB)], Int)
    end

    m =  ModelD(model, x, d, gA, gB)

    loop_cons = handle_loops!(gA, gB, x)
    @objective(m.model, Max, sum(m.d) + loop_cons )


    add_injection(m)
    add_adjacency(m)

    optimize!(m.model)
    display_solution(m, 2)

    return m, 2
end
