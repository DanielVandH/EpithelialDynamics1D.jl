using OrdinaryDiffEq
using ..EpithelialDynamics1D
const CS = EpithelialDynamics1D
using StableRNGs
using Random

## Setting up the problem
s, k, η, t0, t1 = 1.1832, 0.293291, 1.332921, 0.2, 5.0
secs = 15.0
saveat = (t1 - t0) / (24secs)
saveat = t0:saveat:t1
fix_left = true
fix_right = false
a, b = 0.0, 20.0
points = [0.0, 0.2, 0.5, 1.0, 1.2, 1.3, 1.4, 2.0, 2.3, 5.0, 6.7, 8.3, 13.3, 13.4, 14.5, 16.6, 18.8, 19.5, 20.0]
F = (d, p) -> p.k * (p.s - d)
Fp = (k=k, s=s)
G = (d, p) -> p.β * d / p.s
Δt, β = 1e-4, 1e-3
Gp = (β=β, s=s)
prob = CellProblem(;
    initial_time=t0,
    final_time=t1,
    force_law=F,
    force_law_parameters=Fp,
    proliferation_law=G,
    proliferation_law_parameters=Gp,
    initial_condition=points,
    fix_left,
    fix_right,
    damping_constant=η,
    proliferation_period=Δt)
integrator = init(
    prob,
    Tsit5();
    proliferation=true,
    rng=Random.default_rng()
)

## Setting up Gvec 
Gvec = CS.build_Gvec(prob)
@test Gvec == zeros(length(points) - 1)
@inferred CS.build_Gvec(prob)

## Building the proliferation total
Gtot, E = CS.build_proliferation_vec!(prob, Gvec, points)
@inferred CS.build_proliferation_vec!(prob, Gvec, points)
G_manual = zeros(length(points) - 1)
for i in 1:(length(points)-1)
    G_manual[i] = G(abs(points[i+1] - points[i]), Gp)
end
@test E ≈ sum(G_manual) * Δt
@test Gtot ≈ sum(G_manual)
@test Gvec ≈ cumsum(G_manual)
@test Gtot ≈ Gvec[end]

## Testing if a proliferation event occurs 
@test !CS.proliferation_event_occurs(0.5215, 0.3)
@test CS.proliferation_event_occurs(0.6144, 0.7)
@test !CS.proliferation_event_occurs(0.6144, 0.3)
@inferred CS.proliferation_event_occurs(0.5215, 0.3)

## Selecting a cell to proliferate 
_Gvec = [0.1, 0.12, 0.15, 0.18, 0.22, 0.27, 0.5, 0.6, 0.65, 0.71, 0.91, 0.99]
i = CS.select_cell_to_proliferate(0.2, _Gvec)
@test i == 4 + 1
@inferred CS.select_cell_to_proliferate(0.2, _Gvec)
i = CS.select_cell_to_proliferate(0.02, _Gvec)
@test i == 0 + 1
i = CS.select_cell_to_proliferate(0.5, _Gvec)
@test i == 7 + 1
i = CS.select_cell_to_proliferate(0.28, _Gvec)
@test i == 6 + 1

## Test that our method for selecting cells to proliferation is working correctly 
selected = zeros(Int64, length(Gvec))
ns = 100_000
for _ in 1:ns
    local i
    u = Gtot * rand()
    i = CS.select_cell_to_proliferate(u, Gvec)
    selected[i] += 1
end
@test selected / ns ≈ G_manual / Gtot rtol = 1e-1

## Test that we correctly split the cells
orig_pts = deepcopy(points)
CS.split_cell!(orig_pts, 7)
@test length(orig_pts) == length(points) + 1
@test orig_pts[8] == (points[7] + points[8]) / 2
@test issorted(orig_pts)
@test orig_pts[9] == points[8]
@test orig_pts[7] == points[7]
CS.split_cell!(orig_pts, 1)
CS.split_cell!(orig_pts, length(orig_pts) - 1)
@test length(orig_pts) == length(points) + 3
@test orig_pts[2] == (points[1] + points[2]) / 2

## Testing the proliferation affect doesn't break - test it formerly later by comparing to PDE
Gvec = [0.0033 0.0667 0.2039 0.2091 0.2944 0.3389 0.3748 0.4106 0.5064 0.5111 0.5276 0.5892 0.6054 0.7395 0.8092 0.9156 0.9215 0.9363] |> vec
CS.proliferation_affect!(integrator, Gvec)