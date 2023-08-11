function build_proliferation_vec!(prob, Gvec::AbstractVector{T}, nodes) where {T}
    E = zero(T)
    Gtot = zero(T)
    Δt = prob.proliferation_period
    Glaw = prob.proliferation_law
    Gp = prob.proliferation_law_parameters
    for i in eachindex(Gvec) # Gvec is (n-1)×1, nodes is n×1
        Gᵢ = Glaw(nodes[i+1] - nodes[i], Gp)
        Gᵢ = max(Gᵢ, zero(Gᵢ))
        E += Gᵢ * Δt
        Gtot += Gᵢ
        Gvec[i] = Gtot
    end
    return Gtot, E
end

proliferation_event_occurs(u, E) = u < E
select_cell_to_proliferate(u, Gvec) = searchsortedlast(Gvec, u) + 1 # + 1 because we want to avoid i = 0 - note that the last entry will be 1 (or sum(Gvec)), so it's never reached

function split_cell!(nodes, i)
    xᵢ = nodes[i]
    xᵢ₊₁ = nodes[i+1]
    insert!(nodes, i + 1, 0.5(xᵢ + xᵢ₊₁))
    return nothing
end

function proliferation_affect!(integrator, Gvec, rng=Random.default_rng())
    prob = integrator.p
    Gtot, E = build_proliferation_vec!(prob, Gvec, integrator.u)
    u = rand(rng)
    proliferate = proliferation_event_occurs(u, E)
    if !proliferate
        u_modified!(integrator, false)
        return nothing
    else
        u = Gtot * rand(rng)
        i = select_cell_to_proliferate(u, Gvec)
        split_cell!(integrator.u, i)
        resize!(integrator, length(integrator.u))
        resize!(Gvec, length(Gvec) + 1)
        return nothing
    end
    return nothing
end

function build_Gvec(prob::CellProblem)
    n = length(prob.initial_condition)
    Gvec = zeros(eltype(prob.initial_condition), n - 1)
    return Gvec
end

function build_proliferation_callback(prob::CellProblem, rng=Random.default_rng())
    Gvec = build_Gvec(prob)
    f = let Gvec = Gvec, rng = rng
        integrator -> proliferation_affect!(integrator, Gvec, rng)
    end
    Δt = prob.proliferation_period
    callback = PeriodicCallback(f, Δt; save_positions=(false, false))
    return callback
end