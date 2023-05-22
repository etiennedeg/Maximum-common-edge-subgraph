
mutable struct ModelRAD <: AbstractMCESModel
    model::Model
    x::Array{VariableRef,2}
    s::Array{VariableRef,2}
    gA::Graph
    gB::Graph
end

@kwdef struct ConfigRAD <: Config
    name::String
    optimizer::String = "Gurobi"
    max_runtime::Int = 1000
    relaxation::Bool=false
    formulation::String = "RAD"
    adjA::Bool = true
    adjB::Bool = true
end

function add_adjacency_a(m::ModelRAD)
    for u1 in vertices(m.gA), v2 in vertices(m.gB)
        @constraint(m.model,
            sum(m.x[u2, v2] for u2 in neighbors(m.gA, u1))
            - sum(m.x[u1, v1] for v1 in neighbors(m.gB, v2))
            <= m.s[u1, v2]
        )
    end
end

function add_adjacency_b(m::ModelRAD)
    for u1 in vertices(m.gA), v2 in vertices(m.gB)
        @constraint(m.model,
            sum(m.x[u1, v1] for v1 in neighbors(m.gB, v2))
            - sum(m.x[u2, v2] for u2 in neighbors(m.gA, u1))
            <= m.s[u1, v2]
        )
    end
end

function run_config(cf::ConfigRAD, gA::AbstractGraph, gB::AbstractGraph)
    model = get_model(cf)

    # clean_graph!(gA)
    # clean_graph!(gB)
    # gA must be the smallest graph
    if nv(gB) < nv(gA)
        gA, gB = gB, gA
    end

    # if nv(gA) < nv(gB)
    #     add_vertices!(gA, nv(gB)-nv(gA))
    # end

    if cf.relaxation
        @variable(model, 0<=x[i=vertices(gA),j=vertices(gB)]<=1)
        @variable(model, 0<=s[i=vertices(gA),j=vertices(gB)]<=1)
    else
        @variable(model, x[i=vertices(gA),j=vertices(gB)], Bin)
        @variable(model, s[i=vertices(gA),j=vertices(gB)], Bin)
    end

    m = ModelRAD(model, x, s, gA, gB)


    loop_cons = handle_loops!(gA, gB, x)

    startObj = 2*( (cf.adjA ? ne(m.gA) : 0) + (cf.adjB ? ne(m.gB) : 0) )
    obj = startObj - sum(m.s[u, v] for u in vertices(m.gA), v in vertices(m.gB)) + 2*loop_cons
    @objective(m.model, Max, obj)

    add_injection(m)
    cf.adjA && add_adjacency_a(m)
    cf.adjB && add_adjacency_b(m)

    optimize!(m.model)
    coeff = (cf.adjA ? 2 : 0) + (cf.adjB ? 2 : 0)
    display_solution(m, coeff)

    return m, coeff
end
