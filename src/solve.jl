function default_prob_func(prob, proliferation, rng)
    prob_func = if proliferation
        let prob = prob, rng = rng
            (ode_problem, i, repeat) -> begin
                remake(ode_problem; callback=build_proliferation_callback(prob, rng))
            end
        end
    else
        (ode_problem, i, repeat) -> (ode_problem)
    end
    return prob_func
end

function SciMLBase.ODEProblem(prob::CellProblem;
    specialization::Type{S}=SciMLBase.AutoSpecialize,
    jac_prototype=jacobian_sparsity(prob),
    proliferation=false,
    rng=Random.default_rng(),
    sort=false,
    kwargs...) where {S}
    initial_time = prob.initial_time
    final_time = prob.final_time
    time_span = (initial_time, final_time)
    initial_condition = prob.initial_condition
    sort_callback = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> sort!(integrator.u),
        save_positions=(false, false),
        initialize=(c, u, t, integrator) -> c.affect!(integrator),
        finalize=(c, u, t, integrator) -> c.affect!(integrator))
    if proliferation
        f = ODEFunction{true,S}(cell_odes!)
        callback = build_proliferation_callback(prob, rng)
        if sort
            callback = CallbackSet(sort_callback, callback)
        end
        ode_problem = ODEProblem{true,S}(f, initial_condition, time_span, prob; callback=callback, kwargs...)
    else
        f = ODEFunction{true,S}(cell_odes!; jac_prototype=jac_prototype)
        if sort
            ode_problem = ODEProblem{true,S}(f, initial_condition, time_span, prob; callback=sort_callback, kwargs...)
        else
            ode_problem = ODEProblem{true,S}(f, initial_condition, time_span, prob; kwargs...)
        end
    end
    return ode_problem
end
function SciMLBase.NonlinearProblem(prob::SteadyCellProblem; proliferation=false, rng=Random.default_rng(), sort=false, kwargs...)
    proliferation && @warn "Cannot compute the steady state of problems with proliferation. Ignoring the proliferation kwarg."
    ode_prob = ODEProblem(prob.prob; rng, sort, kwargs...)
    nl_prob = NonlinearProblem{true}(ode_prob.f, ode_prob.u0, ode_prob.p; kwargs...)
    return nl_prob
end
function SciMLBase.EnsembleProblem(prob::CellProblem; proliferation=true, rng=Random.default_rng(), sort=false,
    prob_func=default_prob_func(prob, proliferation, rng), kwargs...)
    return EnsembleProblem(ODEProblem(deepcopy(prob); proliferation, rng, sort, kwargs...); prob_func, kwargs...)
end

CommonSolve.init(prob::CellProblem, alg; proliferation=false, rng=Random.default_rng(), sort=false, kwargs...) = CommonSolve.init(ODEProblem(prob; proliferation, rng, sort, kwargs...), alg; kwargs...)
CommonSolve.solve(prob::SteadyCellProblem, alg; proliferation=false, rng=Random.default_rng(), sort=false, kwargs...) = CommonSolve.solve(NonlinearProblem(prob; proliferation, rng, sort, kwargs...), alg; kwargs...)
