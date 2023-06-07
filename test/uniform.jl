using OrdinaryDiffEq
using ..EpithelialDynamics1D
using SteadyStateDiffEq
const CS = EpithelialDynamics1D

force_law = (δ, p) -> p.k * (p.s - δ)
force_law_parameters = (k=1.3, s=1.0)
final_time = 1.0
damping_constant = 1.0
N = 31
initial_condition = (collect ∘ LinRange)(0, 30, N)
prob = CellProblem(;
    force_law,
    force_law_parameters,
    final_time,
    damping_constant,
    initial_condition)
sol = solve(prob, Tsit5(), saveat=0.01)
sprob = SteadyCellProblem(prob)
ssol = solve(sprob, DynamicSS(TRBDF2()))
@test all(≈(ssol.u), sol.u)
@test length(sol) == 101

cell_q = cell_densities.(sol.u)
node_q = node_densities.(sol.u)
@test all(≈(1), reduce(vcat, cell_q))
@test all(≈(1), reduce(vcat, node_q))

cellmpt = map(sol) do u 
    [0.5 * (u[i] +u[i+1]) for i in 1:(length(u)-1)]
end
@test cellmpt ≈ cell_midpoints.(sol.u)