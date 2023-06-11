"""
    continuum_limit(prob::CellProblem,
        mesh_points=copy(prob.initial_condition);
        proliferation=false)

Returns the `FVMProblem` (if the right boundary is fixed) or 
the `MBProblem` (if the right boundary is free) that defines the 
continuum limit to the [`CellProblem`](@ref) `prob`. The argument 
`mesh_points` defines the mesh points to use for the PDE problem, 
which can be a grid of points or a number of points to use. If the problem 
should be considered to have proliferation, set `proliferation = true`.
"""
function continuum_limit(prob::CellProblem,
    mesh_points=copy(prob.initial_condition);
    proliferation=false)
    if prob.fix_right
        return FVMProblem(prob, mesh_points; proliferation)
    else
        return MBProblem(prob, mesh_points; proliferation)
    end
end

"""
    FVMProblem(prob::CellProblem,
        mesh_points=copy(prob.initial_condition);
        diffusion_function=:continuum,
        diffusion_parameters=nothing,
        reaction_function=:continuum,
        reaction_parameters=nothing,
        proliferation=false)

Constructs an `FVMProblem` from a given [`CellProblem`](@ref).
"""
function FiniteVolumeMethod1D.FVMProblem(
    prob::CellProblem,
    mesh_points=copy(prob.initial_condition);
    diffusion_function=:continuum,
    diffusion_parameters=nothing,
    diffusion_theta=nothing,
    reaction_function=:continuum,
    reaction_parameters=nothing,
    reaction_theta=nothing,
    proliferation=false
)
    if !prob.fix_left || prob.initial_condition[begin] ≠ 0.0
        error("The left node must be fixed to zero.")
    end
    if mesh_points isa Integer
        mesh_points = LinRange(prob.initial_condition[begin], prob.initial_condition[end], mesh_points)
    end
    q = _continuum_initial_condition(prob, mesh_points)
    if !isnothing(diffusion_theta)
        diffusion_parameters = (
            θ=diffusion_theta,
            p=diffusion_parameters
        )
    end
    if !isnothing(reaction_theta)
        reaction_parameters = (
            θ=reaction_theta,
            p=reaction_parameters,
        )
    end
    if diffusion_function == :continuum
        D, Dp = _continuum_diffusion_function(prob)
    else
        D, Dp = diffusion_function, diffusion_parameters
    end
    if reaction_function == :continuum
        R, Rp = _continuum_reaction_function(prob, proliferation)
    else
        R, Rp = reaction_function, reaction_parameters
    end
    lhs = FiniteVolumeMethod1D.Neumann(0.0)
    rhs = FiniteVolumeMethod1D.Neumann(0.0)
    return FiniteVolumeMethod1D.FVMProblem(mesh_points, lhs, rhs;
        diffusion_function=D,
        diffusion_parameters=Dp,
        reaction_function=R,
        reaction_parameters=Rp,
        initial_condition=q,
        initial_time=prob.initial_time,
        final_time=prob.final_time)
end

"""
    MBProblem(prob::CellProblem,
        mesh_points=copy(prob.initial_condition);
        diffusion_function=:continuum,
        diffusion_parameters=nothing,
        reaction_function=:continuum,
        reaction_parameters=nothing,
        moving_boundary_function=:continuum,
        moving_boundary_parameters=nothing,
        rhs_function = :continuum,
        rhs_parameters=nothing,
        proliferation=false)

Constructs an `MBProblem` from a given [`CellProblem`](@ref). 
"""
function MovingBoundaryProblems1D.MBProblem(
    prob::CellProblem,
    mesh_points=copy(prob.initial_condition);
    diffusion_function=:continuum,
    diffusion_parameters=nothing,
    diffusion_theta=nothing,
    reaction_function=:continuum,
    reaction_parameters=nothing,
    reaction_theta=nothing,
    moving_boundary_function=:continuum,
    moving_boundary_parameters=nothing,
    moving_boundary_theta=nothing,
    rhs_function=:continuum,
    rhs_parameters=nothing,
    rhs_theta=nothing,
    proliferation=false
)
    if !prob.fix_left || prob.initial_condition[begin] ≠ 0.0
        error("The left node must be fixed to zero.")
    end
    if mesh_points isa Integer
        mesh_points = LinRange(prob.initial_condition[begin], prob.initial_condition[end], mesh_points)
    end
    q = _continuum_initial_condition(prob, mesh_points)
    if !isnothing(diffusion_theta)
        diffusion_parameters = (
            θ=diffusion_theta,
            p=diffusion_parameters
        )
    end
    if !isnothing(reaction_theta)
        reaction_parameters = (
            θ=reaction_theta,
            p=reaction_parameters,
        )
    end
    if !isnothing(moving_boundary_theta)
        moving_boundary_parameters = (
            θ=moving_boundary_theta,
            p=moving_boundary_parameters
        )
    end
    if !isnothing(rhs_theta)
        rhs_parameters = (
            θ=rhs_theta,
            p=rhs_parameters
        )
    end
    if diffusion_function == :continuum
        D, Dp = _continuum_diffusion_function(prob)
    else
        D, Dp = diffusion_function, diffusion_parameters
    end
    if reaction_function == :continuum
        R, Rp = _continuum_reaction_function(prob, proliferation)
    else
        R, Rp = reaction_function, reaction_parameters
    end
    if moving_boundary_function == :continuum
        mb = _continuum_moving_boundary(D, Dp)
    else
        mb = MovingBoundaryProblems1D.Robin(moving_boundary_function, moving_boundary_parameters)
    end
    if rhs_function == :continuum
        rhs = _continuum_rhs(prob, D, Dp)
    else
        rhs = MovingBoundaryProblems1D.Neumann(rhs_function, rhs_parameters)
    end
    lhs = MovingBoundaryProblems1D.Neumann(0.0)
    return MovingBoundaryProblems1D.MBProblem(
        mesh_points ./ mesh_points[end], lhs, rhs, mb;
        diffusion_function=D,
        diffusion_parameters=Dp,
        reaction_function=R,
        reaction_parameters=Rp,
        initial_condition=q,
        initial_time=prob.initial_time,
        initial_endpoint=prob.initial_condition[end],
        final_time=prob.final_time
    )
end

function _continuum_initial_condition(prob, mesh_points)
    q = node_densities(prob.initial_condition)
    interp = LinearInterpolation{true}(q, prob.initial_condition)
    pde_initial_condition = interp.(mesh_points)
    return pde_initial_condition
end

function _continuum_diffusion_function(prob)
    F′ = let cell_prob = prob
        q -> begin
            F = cell_prob.force_law
            Fp = cell_prob.force_law_parameters
            ForwardDiff.derivative(s -> F(s, Fp), inv(q))
        end
    end
    Dp = (η⁻¹=inv(prob.damping_constant), F′=F′)
    D = (q, x, t, p) -> -p.η⁻¹ / q^2 * p.F′(q)
    return D, Dp
end

function _continuum_reaction_function(prob, proliferation)
    if proliferation
        G = prob.proliferation_law
        Gp = prob.proliferation_law_parameters
        Rp = (G=G, Gp=Gp)
        R = (q, x, t, p) -> q * p.G(inv(q), p.Gp)
        return R, Rp
    else
        Rp = nothing
        R = (q, x, t, p) -> zero(q)
        return R, Rp
    end
end

function _continuum_rhs(prob, diffusion, diffusion_parameters)
    η = prob.damping_constant
    η⁻¹ = inv(η)
    F = prob.force_law
    Fp = prob.force_law_parameters
    rhs_f = (q, t, p) -> -2q * p.η⁻¹ * F(inv(q), Fp) * inv(p.D(q, nothing, t, p.Dp))
    rhs_p = (η⁻¹=η⁻¹, F=F, Fp=Fp, D=diffusion, Dp=diffusion_parameters)
    return MovingBoundaryProblems1D.Neumann(rhs_f, rhs_p)
end

function _continuum_moving_boundary(diffusion, diffusion_parameters)
    mb_f = (q, t, p) -> (zero(q), -p.D(q, nothing, t, p.Dp) * inv(q))
    mb_p = (D=diffusion, Dp=diffusion_parameters)
    return MovingBoundaryProblems1D.Robin(mb_f, mb_p)
end

