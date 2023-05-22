
mutable struct ModelML <: AbstractMCESModel
    model::Model
    x::Array{VariableRef,2}
    y::Array{VariableRef,2}
    gA::Graph
    gB::Graph
end

@kwdef struct ConfigML <: Config
    name::String
    optimizer::String="Gurobi"
    max_runtime::Int = 1000
    relaxation::Bool=false
    formulation::String = "M"
    degreeCuts::Bool=false
end

function add_y_bounds(m::ModelML)
    for u1 in vertices(m.gA)
        for u2 in vertices(m.gA)
            if (u2 < u1) || !has_edge(m.gA, u1, u2)
                fix(m.y[u1, u2], 0; force = true)
            end
        end
    end
end

function add_adjacency(m::ModelML)
    for u1 in vertices(m.gA), u2 in outneighbors(m.gA, u1)
        for v1 in vertices(m.gB)
            @constraint(m.model, m.y[min(u1,u2), max(u1,u2)] + m.x[u1,v1] <= 1 + sum(m.x[u2,v2] for v2 in outneighbors(m.gB, v1)) )
        end
    end
end

function add_degree_cuts(m::ModelML)
    for u1 in vertices(m.gA)
        @constraint(m.model, sum(m.y[min(u1, u2), max(u1, u2)] for u2 in neighbors(m.gA, u1)) <= sum(min(degree(m.gA, u1), degree(m.gB, v1))*m.x[u1, v1] for v1 in vertices(m.gB)) )
    end
end

function run_config(cf::ConfigML, gA::AbstractGraph, gB::AbstractGraph)
    model = get_model(cf)

    clean_graph!(gA)
    clean_graph!(gB)

    # gA must be the smallest graph
    if nv(gB) < nv(gA)
        gA, gB = gB, gA
    end

    if cf.relaxation
        @variable(model, 0 <= x[i=vertices(gA),j=vertices(gB)] <= 1)
        @variable(model, 0 <= y[i=vertices(gA),j=vertices(gA)] <= 1)
    else
        @variable(model, x[i=vertices(gA),j=vertices(gB)], Bin)
        @variable(model, y[i=vertices(gA),j=vertices(gA)], Bin)
    end

    m = ModelML(model, x, y, gA, gB)

    loop_cons = handle_loops!(gA, gB, x)
    @objective(m.model, Max, sum(m.y) + loop_cons)

    add_y_bounds(m)
    add_injection(m)
    add_adjacency(m)
    cf.degreeCuts && add_degree_cuts(m)

    optimize!(m.model)
    display_solution(m, 1)

    return m, 1
end
