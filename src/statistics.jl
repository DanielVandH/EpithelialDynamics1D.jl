"""
    cell_densities(cell_positions::AbstractVector{T}) where {T<:Number}

Compute the cell densities from the cell positions. The returned vector `q` has 
`length(q) = length(cell_positions) - 1`, with `q[i]` the density of the
`i`th cell `(cell_positions[i], cell_positions[i+1])` given by 
`1/(cell_positions[i+1] - cell_positions[i])`.

If you want estimates of the densities at each node rather than for a cell, 
see [`node_densities`](@ref).
"""
function cell_densities(cell_positions::AbstractVector{T}) where {T<:Number}
    q = similar(cell_positions, length(cell_positions) - 1)
    for i in eachindex(q)
        q[i] = 1 / (cell_positions[i+1] - cell_positions[i])
    end
    return q
end

"""
    cell_midpoints(cell_positions::AbstractVector{T}) where {T<:Number}

Compute the midpoints of the cells from the cell positions. The returned vector
`x` has `length(x) = length(cell_positions) - 1`, with `x[i]` the midpoint of the
`i`th cell `(cell_positions[i], cell_positions[i+1])` given by 
`0.5(cell_positions[i] + cell_positions[i+1])`.
"""
function cell_midpoints(cell_positions::AbstractVector{T}) where {T<:Number}
    x = similar(cell_positions, length(cell_positions) - 1)
    for i in eachindex(x)
        x[i] = 0.5(cell_positions[i] + cell_positions[i+1])
    end
    return x
end

"""
    node_densities(cell_positions::AbstractVector{T}) where {T<:Number}

Compute the cell densities from the cell positions, assigning a density to each node. 
The `i`th density is given by `2/(cell_positions[i+1] - cell_positions[i-1])` if
`i` is not the first or last node, `1/(cell_positions[i+1] - cell_positions[i])`
if `i` is the first node, and `1/(cell_positions[i] - cell_positions[i-1])` if `i`
is the last node.

If you want estimates of the densities for each cell rather than at each node,
see [`cell_densities`](@ref).
"""
function node_densities(cell_positions::AbstractVector{T}) where {T<:Number}
    q = similar(cell_positions)
    for i in eachindex(q)
        if i == firstindex(q)
            q[i] = 1 / (cell_positions[i+1] - cell_positions[i])
        elseif i == lastindex(q)
            q[i] = 1 / (cell_positions[i] - cell_positions[i-1])
        else
            q[i] = 2 / (cell_positions[i+1] - cell_positions[i-1])
        end
    end
    return q
end

"""
    get_knots(sol, num_knots = 500; indices = eachindex(sol), use_max=true)

Computes knots for each time, covering the extremum of the cell positions across all 
cell simulations. You can restrict the simultaions to consider using the `indices`.
If `use_max` is true, then the knots will be obtained by taking the extreme node positions 
for each `time`, otherwise the average is used.
"""
function get_knots(sol::EnsembleSolution, num_knots=500; indices=eachindex(sol), use_extrema=true)
    @static if VERSION < v"1.7"
        knots = Vector{LinRange{Float64}}(undef, length(first(sol)))
    else
        knots = Vector{LinRange{Float64,Int}}(undef, length(first(sol)))
    end
    times = first(sol).t
    Base.Threads.@threads for i in eachindex(times)
        local a, b
        if use_extrema
            a = Inf
            b = -Inf
        else
            a = 0.0
            b = 0.0
            ctr = 0
        end
        for j in indices
            _a = sol[j][i][begin]
            _b = sol[j][i][end]
            if use_extrema
                a = min(a, _a)
                b = max(b, _b)
            else
                a += _a
                b += _b
                ctr += 1
            end
        end
        if !use_extrema
            a /= ctr
            b /= ctr
        end
        knots[i] = LinRange(a, b, num_knots)
    end
    return knots
end
function get_knots(sol::ODESolution, num_knots=500)
    knots = map(sol) do r
        LinRange(r[begin], r[end], num_knots)
    end
    return knots
end

"""
    node_densities(sol::EnsembleSolution; num_knots=500, knots=get_knots(sol, num_knots), alpha=0.05, interp_fnc=(u, t) -> LinearInterpolation{true}(u, t))

Computes summary statistics for the node densities from an `EnsembleSolution` to a [`CellProblem`](@ref).

# Arguments 
- `sol::EnsembleSolution`: The ensemble solution to a `CellProblem`.

# Keyword Arguments
- `indices = eachindex(sol)`: The indices of the cell simulations to consider. 
- `num_knots::Int = 500`: The number of knots to use for the spline interpolation.
- `use_extrema::Bool = true`: Whether to use the extrema of the cell positions for the knots, or the average.
- `knots::Vector{Vector{Float64}} = get_knots(sol, num_knots; indices, use_extrema)`: The knots to use for the spline interpolation.
- `alpha::Float64 = 0.05`: The significance level for the confidence intervals.
- `interp_fnc = (u, t) -> LinearInterpolation{true}(u, t)`: The function to use for constructing the interpolant.

# Outputs 
- `q::Vector{Vector{Vector{Float64}}}`: The node densities for each cell simulation.
- `r::Vector{Vector{Vector{Float64}}}`: The cell positions for each cell simulation.
- `means::Vector{Vector{Float64}}`: The mean node densities for each cell simulation.
- `lowers::Vector{Vector{Float64}}`: The lower bounds of the confidence intervals for the node densities for each cell simulation.
- `uppers::Vector{Vector{Float64}}`: The upper bounds of the confidence intervals for the node densities for each cell simulation.
- `knots::Vector{Vector{Float64}}`: The knots used for the spline interpolation.
"""
function node_densities(sol::EnsembleSolution;
    indices=eachindex(sol),
    num_knots=500,
    use_extrema=true,
    knots=get_knots(sol, num_knots; indices, use_extrema),
    alpha=0.05,
    interp_fnc=(u, t) -> LinearInterpolation{true}(u, t))
    q = Vector{Vector{Vector{Float64}}}(undef, length(indices))
    r = Vector{Vector{Vector{Float64}}}(undef, length(indices))
    Base.Threads.@threads for i in eachindex(indices)
        q[i] = node_densities.(sol[indices[i]].u)
        r[i] = sol[indices[i]].u
    end
    nt = length(first(sol))
    nsims = length(indices)
    q_splines = zeros(num_knots, nt, nsims)
    q_means = [zeros(num_knots) for _ in 1:nt]
    q_lowers = [zeros(num_knots) for _ in 1:nt]
    q_uppers = [zeros(num_knots) for _ in 1:nt]
    Base.Threads.@threads for k in 1:nsims
        for j in 1:nt
            densities = q[k][j]
            cell_positions = r[k][j]
            interp = interp_fnc(densities, cell_positions)
            for i in eachindex(knots[j])
                if knots[j][i] > r[k][j][end]
                    q_splines[i, j, k] = 0.0
                else
                    q_splines[i, j, k] = max(0.0, interp(knots[j][i]))
                end
            end
        end
    end
    Base.Threads.@threads for j in 1:nt
        knot_range = knots[j]
        for i in eachindex(knot_range)
            q_values = @views q_splines[i, j, :]
            q_means[j][i] = mean(q_values)
            q_lowers[j][i] = quantile(q_values, alpha / 2)
            q_uppers[j][i] = quantile(q_values, 1 - alpha / 2)
        end
    end
    return (q=q, r=r, means=q_means, lowers=q_lowers, uppers=q_uppers, knots=knots)
end

"""
    cell_numbers(sol::EnsembleSolution; indices = eachindex(sol), alpha=0.05)

Computes summary statistics for the cell numbers from an `EnsembleSolution` to a [`CellProblem`](@ref).

# Arguments
- `sol::EnsembleSolution`: The ensemble solution to a `CellProblem`.

# Keyword Arguments
- `indices = eachindex(sol)`: The indices of the cell simulations to consider.
- `alpha::Float64 = 0.05`: The significance level for the confidence intervals.

# Outputs
- `N::Vector{Vector{Int}}`: The cell numbers for each cell simulation.
- `means::Vector{Float64}`: The mean cell numbers for each cell simulation.
- `lowers::Vector{Float64}`: The lower bounds of the confidence intervals for the cell numbers for each cell simulation.
- `uppers::Vector{Float64}`: The upper bounds of the confidence intervals for the cell numbers for each cell simulation.
"""
function cell_numbers(sol::EnsembleSolution; indices=eachindex(sol), alpha=0.05)
    N = map(indices) do i
        length.(sol[i].u) .- 1
    end |> x -> reduce(hcat, x)
    N_means = zeros(size(N, 1))
    N_lowers = zeros(size(N, 1))
    N_uppers = zeros(size(N, 1))
    for (i, n) in (enumerate ∘ eachrow)(N)
        N_means[i] = mean(n)
        N_lowers[i] = quantile(n, alpha / 2)
        N_uppers[i] = quantile(n, 1 - alpha / 2)
    end
    @static if VERSION ≥ v"1.9"
        return (N=convert(Vector{Vector{Int}}, eachcol(N)), means=N_means, lowers=N_lowers, uppers=N_uppers)
    else
        return (N=convert(Vector{Vector{Int}}, (collect ∘ eachcol)(N)), means=N_means, lowers=N_lowers, uppers=N_uppers)
    end
end

"""
    leading_edges(sol::EnsembleSolution; indices = eachindex(sol), alpha=0.05)

Computes summary statistics for the leading edges from an `EnsembleSolution` to a [`CellProblem`](@ref).

# Arguments
- `sol::EnsembleSolution`: The ensemble solution to a `CellProblem`.

# Keyword Arguments
- `indices = eachindex(sol)`: The indices of the cell simulations to consider.
- `alpha::Float64 = 0.05`: The significance level for the confidence intervals.

# Outputs
- `L::Vector{Vector{Float64}}`: The leading edges for each cell simulation.
- `means::Vector{Float64}`: The mean leading edges for each cell simulation.
- `lowers::Vector{Float64}`: The lower bounds of the confidence intervals for the leading edges for each cell simulation.
- `uppers::Vector{Float64}`: The upper bounds of the confidence intervals for the leading edges for each cell simulation.
"""
function leading_edges(sol::EnsembleSolution; indices=eachindex(sol), alpha=0.05)
    L = map(indices) do i
        map(sol[i]) do _sol
            _sol[end]
        end
    end |> x -> reduce(hcat, x)
    L_means = zeros(size(L, 1))
    L_lowers = zeros(size(L, 1))
    L_uppers = zeros(size(L, 1))
    for (i, n) in (enumerate ∘ eachrow)(L)
        L_means[i] = mean(n)
        L_lowers[i] = quantile(n, alpha / 2)
        L_uppers[i] = quantile(n, 1 - alpha / 2)
    end
    @static if VERSION ≥ v"1.9"
        return (L=convert(Vector{Vector{Float64}}, eachcol(L)), means=L_means, lowers=L_lowers, uppers=L_uppers)
    else
        return (L=convert(Vector{Vector{Float64}}, (collect ∘ eachcol)(L)), means=L_means, lowers=L_lowers, uppers=L_uppers)
    end
end

"""
    integrate_pde(sol, i)

Integrate the PDE solution `sol` at time `i`, giving an approximation of the number of cells at `sol.t[i]`.
"""
function integrate_pde(sol, i)
    prob = sol.prob.p
    x = prob.geometry.mesh_points
    u = sol.u[i]
    if typeof(prob) <: FVMProblem
        interp = LinearInterpolation{true}(u, x)
        N = DataInterpolations.integral(interp, x[1], x[end])
        return N
    else # if typeof(prob) <: MBProblem 
        x = prob.geometry.mesh_points
        interp = LinearInterpolation{true}(@views(u[begin:(end-1)]), x)
        N = DataInterpolations.integral(interp, zero(eltype(x)), one(eltype(x)))
        L = u[end]
        return L * N
    end
end
