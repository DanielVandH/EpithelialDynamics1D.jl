using ..EpithelialDynamics1D

force_law = (δ, p) -> p.k * (p.s - δ)
force_law_parameters = (k=1.3, s=0.7)
proliferation_law = (δ, p) -> p.β
proliferation_law_parameters = (β = 0.1,)
proliferation_period = 1e-2
initial_time = 0.0
final_time = 1.0
fix_left = true
fix_right = false
damping_constant = 1.2
initial_condition = [0.0, 0.5, 1.0]
prob = CellProblem(;
    force_law=force_law,
    force_law_parameters=force_law_parameters,
    proliferation_law=proliferation_law,
    proliferation_law_parameters=proliferation_law_parameters,
    proliferation_period=proliferation_period,
    initial_time=initial_time,
    final_time=final_time,
    fix_left=fix_left,
    fix_right=fix_right,
    damping_constant=damping_constant,
    initial_condition=initial_condition
)
@test prob.force_law_parameters.k == 1.3
@test prob.force_law_parameters.s == 0.7
@test prob.proliferation_law_parameters.β == 0.1
@test prob.proliferation_period == 1e-2
@test prob.initial_time == 0.0
@test prob.final_time == 1.0
@test prob.fix_left
@test !prob.fix_right
@test prob.damping_constant == 1.2
@test prob.initial_condition == [0.0, 0.5, 1.0]
@test prob.force_law(0.1, prob.force_law_parameters) == 1.3 * (0.7 - 0.1)
@test prob.proliferation_law(0.1, prob.proliferation_law_parameters) == 0.1

sprob = SteadyCellProblem(prob)
@test sprob.prob == prob