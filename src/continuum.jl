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
    if !prob.fix_left || prob.initial_condition[begin] ≠ 0.0
        error("The left node must be fixed to zero.")
    end
    if mesh_points isa Integer
        mesh_points = LinRange(prob.initial_condition[begin], prob.initial_condition[end], mesh_points)
    end
    q = _continuum_initial_condition(prob, mesh_points)
    if prob.fix_right
        return _fvm_continuum_limit(prob, mesh_points, q, proliferation)
    else
        return _mb_continuum_limit(prob, mesh_points, q, proliferation)
    end
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

function _fvm_continuum_limit(prob::CellProblem, mesh_points, q, proliferation)
    D, Dp = _continuum_diffusion_function(prob)
    R, Rp = _continuum_reaction_function(prob,  proliferation)
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

function _mb_continuum_limit(prob::CellProblem, mesh_points, q, proliferation)
    D, Dp = _continuum_diffusion_function(prob)
    R, Rp = _continuum_reaction_function(prob, proliferation)
    lhs = MovingBoundaryProblems1D.Neumann(0.0)
    rhs = _continuum_rhs(prob, D, Dp)
    mb = _continuum_moving_boundary(D, Dp)
    return MovingBoundaryProblems1D.MBProblem(mesh_points ./ mesh_points[end], lhs, rhs, mb;
        diffusion_function=D,
        diffusion_parameters=Dp,
        reaction_function=R,
        reaction_parameters=Rp,
        initial_condition=q,
        initial_time=prob.initial_time,
        initial_endpoint=prob.initial_condition[end],
        final_time=prob.final_time)
end