
mutable struct ModelBMPS <: AbstractMCESModel
    model::Model
    x::Array{VariableRef,2}
    y::JuMP.Containers.DenseAxisArray{VariableRef, 1, Tuple{Vector{Tuple{Int64, Int64, Int64, Int64}}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64, Int64, Int64}, Int64}}}}
    gA::Graph
    gB::Graph
end

Base.@kwdef struct ConfigBMPS <: Config
    name::String
    optimizer::String = "Gurobi"
    max_runtime::Int = 1000
    relaxation::Bool = false
    formulation::String = "BMPS"
    symetrized::Bool = false
end
#
function get_d(m::ModelBMPS, u1, u2, v1, v2)
    minu = min(u1,u2)
    maxu = max(u1,u2)
    minv = min(v1,v2)
    maxv = max(v1,v2)
    return m.y[(minu, maxu, minv, maxv)]
end

function get_d2(m::ModelBMPS, u1, u2, v1, v2)
    if u2 < u1
        return m.y[(u2, u1, v2, v1)]
    else
        return m.y[(u1, u2, v1, v2)]
    end
end

function add_adjacency(m::ModelBMPS)
    for u1 in vertices(m.gA), u2 in outneighbors(m.gA, u1)
        (u1 <= u2) && continue
        for v1 in vertices(m.gB)
            @constraint(m.model, sum(get_d(m, u1,u2,v1,v2) for v2 in outneighbors(m.gB, v1)) <= m.x[u1,v1] + m.x[u2,v1])
        end
    end
    for v1 in vertices(m.gB), v2 in outneighbors(m.gB, v1)
        (v1 <= v2) && continue
        for u1 in vertices(m.gA)
            @constraint(m.model, sum(get_d(m, u1,u2,v1,v2) for u2 in outneighbors(m.gA, u1)) <= m.x[u1,v1] + m.x[u1,v2])
        end
    end
end

function add_adjacency_s(m::ModelBMPS)
    for e1 in edges(m.gA), v1 in vertices(m.gB)
        @constraint(m.model, sum(get_d2(m, e1.src,e1.dst,v1,v2) for v2 in neighbors(m.gB, v1)) <= m.x[e1.src,v1])
        @constraint(m.model, sum(get_d2(m, e1.src,e1.dst,v2,v1) for v2 in neighbors(m.gB, v1)) <= m.x[e1.dst,v1])
    end
    for e2 in edges(m.gB), u1 in vertices(m.gA)
        @constraint(m.model, sum(get_d2(m, u1, u2, e2.src,e2.dst) for u2 in neighbors(m.gA, u1)) <= m.x[u1, e2.src])
        @constraint(m.model, sum(get_d2(m, u2, u1, e2.src,e2.dst) for u2 in neighbors(m.gA, u1)) <= m.x[u1, e2.dst])
    end
end

function run_config(cf::ConfigBMPS, gA::AbstractGraph, gB::AbstractGraph)
    model = get_model(cf)

    clean_graph!(gA)
    clean_graph!(gB)

    # gA must be the smallest graph
    if nv(gB) < nv(gA)
        gA, gB = gB, gA
    end

    if cf.symetrized
        y_indices = [(src(ea), dst(ea), src(eb), dst(eb)) for ea in edges(gA), eb in all_edges(gB) if src(ea) != dst(ea) && src(eb) != dst(eb)]
    else
        y_indices = [(src(ea), dst(ea), src(eb), dst(eb)) for ea in edges(gA), eb in all_edges(gB) if src(ea) != dst(ea) && src(eb) != dst(eb)]
    end
    if cf.relaxation
        @variable(model, 0<=x[i=vertices(gA),j=vertices(gB)]<=1)
        @variable(model, 0<=y[y_indices])
    else
        @variable(model, x[i=vertices(gA),j=vertices(gB)], Bin)
        @variable(model, y[y_indices], Bin)
    end

    m = ModelBMPS(model,x, y, gA, gB)

    loop_cons = handle_loops!(gA, gB, x)
    if !cf.symetrized
        @objective(m.model, Max, sum(m.y[(e1.src, e1.dst, e2.src, e2.dst)] for e1 in edges(m.gA), e2 in edges(m.gB)) + loop_cons)
    else
        @objective(m.model, Max, sum(m.y[(e1.src, e1.dst, e2.src, e2.dst)] + m.y[(e1.src, e1.dst, e2.dst, e2.src)] for e1 in edges(m.gA), e2 in edges(m.gB)) + loop_cons)
    end

    add_injection(m)
    if !cf.symetrized
        add_adjacency(m)
    else
        add_adjacency_s(m)
    end


    optimize!(m.model)
    display_solution(m, 1)

    return m, 1
end
