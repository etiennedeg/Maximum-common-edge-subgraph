
macro timeout(seconds, expr)
    quote
        tsk = @task $expr
        schedule(tsk)
        println(tsk)
        println($seconds)
        Timer($seconds) do timer
            istaskdone(tsk) || Base.throwto(tsk, InterruptException())
        end
        fetch(tsk)
    end
end

"""
`run` launch the solve process in a separate thread, to allow the interruption of the solve
"""
function run_config_in_thread(cf, gA, gB)
    local results
    try
        tsk = @task begin results = run_config(cf, gA, gB) end
        schedule(tsk)
        Timer(cf.max_runtime) do timer
            istaskdone(tsk) || Base.throwto(tsk, TaskFailedException(tsk))
        end
        fetch(tsk)
    catch e
        if isa(e, TaskFailedException)
            if isa(e.task.result, TaskFailedException)
                results = :timeout
            elseif isa(e.task.result, InterruptException)
                rethrow(e.task.result)
            else
                rethrow(e.task.result)
            end
        else
            rethrow(e)
        end
    end
    return results
end


function handle_loops!(gA, gB, x)
	cons = GenericAffExpr{Float64,VariableRef}()
	for u in vertices(gA)
		!has_edge(gA, u, u) && continue
		for v in vertices(gB)
			!has_edge(gB, v, v) && continue
			add_to_expression!(cons, x[u, v])
		end
	end
	free_loops!(gA)
	free_loops!(gB)
	return cons
end

function add_injection(m::AbstractMCESModel)
    @constraint(m.model, [u in vertices(m.gA)], sum(m.x[u, :]) == 1)
    @constraint(m.model, [v in vertices(m.gB)], sum(m.x[:, v]) <= 1)
end

function get_model(cf::Config)
	if cf.optimizer == "Gurobi"
		model = Model(optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => cf.max_runtime))
	elseif cf.optimizer == "CPLEX"
		model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_TILIM" => cf.max_runtime))
	elseif cf.optimizer == "SCIP"
		model = Model(optimizer_with_attributes(SCIP.Optimizer, "limits/time" => cf.max_runtime))
	end
	return model
end

function display_solution(m::AbstractMCESModel, coeff)
	println(termination_status(m.model))
	if termination_status(m.model) == MOI.OPTIMAL
		println("objective value: $(objective_value(m.model)/coeff)")
		for i in vertices(m.gA), j in vertices(m.gB)
			if value(m.x[i,j]) > 0.0000001
				println("x_{$i, $j} = ", value(m.x[i,j]) )
			end
		end
	end
	# if termination_status(m.model) == MOI.OPTIMAL
	# 	println("objective value: $(objective_value(m.model)/coeff)")
	# 	for i in vertices(m.gA), j in vertices(m.gA)
	# 		if value(m.z[i,j]) > 0.0000001
	# 			println("z_{$i, $j} = ", value(m.z[i,j]) )
	# 		end
	# 	end
	# end
end
