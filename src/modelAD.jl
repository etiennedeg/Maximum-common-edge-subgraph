
mutable struct ModelAD <: AbstractMCESModel
    model::Model
    x::Array{VariableRef,2}
    s::Array{VariableRef,2}
    t::Array{VariableRef,2}
    gA::Graph
    gB::Graph
end

@kwdef struct ConfigAD <: Config
    name::String
    optimizer::String = "Gurobi"
    max_runtime::Int = 1000
    relaxation::Bool=false
    formulation::String = "AD"
end

function add_adjacency(m::ModelAD)
    for u2 in vertices(m.gA), v1 in vertices(m.gB)
        @constraint(m.model,
                - sum(m.x[u1, v1] for u1 in neighbors(m.gA, u2))
                + sum(m.x[u2, v2] for v2 in neighbors(m.gB, v1))
                + m.t[u2, v1] - m.s[u2, v1] == 0
        )
    end
end

function run_config(cf::ConfigAD, gA::AbstractGraph, gB::AbstractGraph)
    model = get_model(cf)

    # gA must be the smallest graph
    if nv(gB) < nv(gA)
        gA, gB = gB, gA
    end

    # graphs must be the same size
    if nv(gA) < nv(gB)
        add_vertices!(gA, nv(gB)-nv(gA))
    end

    if cf.relaxation
        @variable(model, 0<=x[i=vertices(gA),j=vertices(gB)]<=1)
        @variable(model, 0<=s[i=vertices(gA),j=vertices(gB)]<=1)
        @variable(model, 0<=t[i=vertices(gA),j=vertices(gB)]<=1)
    else
        @variable(model, x[i=vertices(gA),j=vertices(gB)], Bin)
        @variable(model, s[i=vertices(gA),j=vertices(gB)], Bin)
        @variable(model, t[i=vertices(gA),j=vertices(gB)], Bin)
    end

    m = ModelAD(model, x, s, t, gA, gB)

    loop_cons = handle_loops!(gA, gB, x)
    obj = 2*(ne(m.gA) + ne(m.gB)) - sum(m.s[u, v] + m.t[u, v] for u in vertices(m.gA), v in vertices(m.gB)) + 2*loop_cons
    @objective(m.model, Max, obj)

    add_injection(m)
    add_adjacency(m)

    optimize!(m.model)
    display_solution(m, 4)

    return m, 4
end
