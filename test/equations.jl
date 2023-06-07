using OrdinaryDiffEq
using LinearAlgebra
using SparseArrays
using ..EpithelialDynamics1D
const CS = EpithelialDynamics1D

## Test the cell ODEs
function _cell_odes!(dr::AbstractVector{T}, r, p, t) where {T}
    (; α, s, fix_left, fix_right) = p
    if !fix_left
        dr[1] = α * (r[2] - r[1] - s)
    else
        dr[1] = zero(T)
    end
    if !fix_right
        dr[end] = α * (r[end-1] - r[end] + s)
    else
        dr[end] = zero(T)
    end
    for i in (firstindex(dr)+1):(lastindex(dr)-1)
        dr[i] = α * (r[i-1] - 2r[i] + r[i+1])
    end
    return nothing
end

s, k, η, t0, t1 = 0.3, 0.293291, 1.332921, 0.2, 5.0
secs = 15.0
saveat = (t1 - t0) / (24secs)
saveat = t0:saveat:t1
a, b = 0.0, 20.0
points = LinRange(a, b, 1000) |> collect
for (fix_left, fix_right) in Iterators.product((false, true), (false, true))
    prob = CellProblem(; initial_time=t0, final_time=t1, fix_left, fix_right, initial_condition=points,
        damping_constant=η, force_law=(δ, p) -> p.k * (p.s - δ), force_law_parameters=(k=k, s=s))
    sol = solve(prob, Tsit5(), saveat=saveat)
    for i in eachindex(sol.t)
        q = zero(points)
        dr = zero(q)
        dr2 = zero(q)
        CS.cell_odes!(dr, sol.u[i], prob, sol.t[i])
        _cell_odes!(dr2, sol.u[i], (α=k / η, s=s, fix_left=prob.fix_left, fix_right=prob.fix_right), sol.t[i])
        @test dr ≈ dr2 atol = 1e-6
    end
end

## Test the Jacobian 
prob = CellProblem(; initial_time=t0, final_time=t1, fix_left=true, fix_right=false, initial_condition=points,
damping_constant=η, force_law=(δ, p) -> p.k * (p.s - δ), force_law_parameters=(k=k, s=s))
jac = CS.jacobian_sparsity(prob)
_jac = zeros(size(jac))
for i in eachindex(points)
    itr = if i == firstindex(points) 
        [1, 2]
    elseif i == lastindex(points)
        [length(points)-1,length(points)]
    else
        [i,i-1,i+1]
    end
    _jac[i, itr] .= 1
end
@test jac == _jac 
@test nnz(jac) == sum(_jac) == 3length(points) - 2 
@test Tridiagonal(jac) == jac

