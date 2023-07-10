using OrdinaryDiffEq
using ..EpithelialDynamics1D
using SteadyStateDiffEq
using StatsBase
using FiniteVolumeMethod1D
using LinearSolve
using LinearAlgebra
using CairoMakie
using ReferenceTests
using DataInterpolations
using MovingBoundaryProblems1D
const CS = EpithelialDynamics1D

@static if VERSION < v"1.7"
    stack(x) = reduce(hcat, x)
end

initial_condition = [LinRange(0, 15, 16); LinRange(15, 30, 32)] |> unique!
fig = Figure(fontsize=33)
ax = Axis(fig[1, 1], xlabel=L"x", width=600, height=200)
scatter!(ax, initial_condition, zero(initial_condition), color=:black, markersize=13)
hideydecorations!(ax)
resize_to_layout!(fig)
fig
fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
@test_reference joinpath(fig_path, "step_function_initial_condition.png") fig

@testset "Fixed Boundary" begin
    println("Starting Fixed Boundary example.")
    force_law = (δ, p) -> p.k * (p.s - δ)
    force_law_parameters = (k=10.0, s=0.2)
    final_time = 100.0
    damping_constant = 1.0
    initial_condition = [LinRange(0, 15, 16); LinRange(15, 30, 32)] |> unique!
    prob = CellProblem(;
        force_law,
        force_law_parameters,
        final_time,
        damping_constant,
        initial_condition)
    for sol in (solve(prob, Tsit5(), saveat=10.0), solve(prob, Tsit5(), saveat=10.0, sort=true))
        sol = solve(prob, Tsit5(), saveat=10.0)
        @test all(≈(LinRange(0, 30, 500)), get_knots(sol))

        fvm_prob = continuum_limit(prob, 1000)
        fvm_sol = solve(fvm_prob, TRBDF2(linsolve=KLUFactorization()), saveat=10.0)

        # Test the cell densities 
        cell_densities_sol = cell_densities.(sol.u)
        node_densities_sol = node_densities.(sol.u)
        for i in eachindex(sol)
            r = sol.u[i]
            q = cell_densities_sol[i]
            nq = node_densities_sol[i]
            for j in eachindex(r)
                if j == firstindex(r)
                    @test nq[j] ≈ 1 / (r[j+1] - r[j])
                    @test q[j] ≈ 1 / (r[j+1] - r[j])
                elseif j == lastindex(r)
                    @test nq[j] ≈ 1 / (r[j] - r[j-1])
                else
                    @test nq[j] ≈ 2 / (r[j+1] - r[j-1])
                    @test q[j] ≈ 1 / (r[j+1] - r[j])
                end
            end
        end

        # Compare with the continuum limit 
        let x = fvm_prob.geometry.mesh_points
            fig = Figure(fontsize=36)
            colors = (:red, :blue, :darkgreen, :black, :orange, :magenta, :cyan, :yellow, :brown, :gray, :lightblue)
            ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"q(x, t)",
                width=600, height=300)
            [lines!(ax, sol.u[i], node_densities_sol[i], color=colors[i], linewidth=2) for i in eachindex(sol)]
            [lines!(ax, x, fvm_sol.u[i], color=colors[i], linewidth=4, linestyle=:dashdot) for i in eachindex(sol)]
            resize_to_layout!(fig)
            fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
            @test_reference joinpath(fig_path, "step_function.png") fig
            fig
        end

        # Look at the steady state 
        sprob = SteadyCellProblem(prob)
        ssol = solve(sprob, DynamicSS(TRBDF2()))
        m = mean(diff(ssol.u))
        @test all(x -> isapprox(x, m, rtol=1e-3), diff(ssol.u))
        q = cell_densities(ssol.u)
        @test all(x -> isapprox(x, 1.5333154004, rtol=1e-3), q)
    end
end

@testset "Moving Boundary" begin
    println("Starting Moving Boundary example.")
    force_law = (δ, p) -> p.k * (p.s - δ)
    force_law_parameters = (k=10.0, s=0.2)
    final_time = 500.0
    damping_constant = 1.0
    initial_condition = [LinRange(0, 15, 16); LinRange(15, 30, 32)] |> unique!
    prob = CellProblem(;
        force_law,
        force_law_parameters,
        final_time,
        damping_constant,
        initial_condition,
        fix_right=false)

    for sol in (solve(prob, Tsit5(), saveat=50.0, sort=true), solve(prob, Tsit5(), saveat=50.0))
        sol = solve(prob, Tsit5(), saveat=50.0)
        _knots = get_knots(sol, 700)
        for i in eachindex(sol)
            @test _knots[i] ≈ LinRange(sol.u[i][begin], sol.u[i][end], 700)
        end

        mb_prob = continuum_limit(prob, 1000)
        mb_sol = solve(mb_prob, TRBDF2(linsolve=KLUFactorization()), saveat=50.0)

        # Test the cell densities 
        cell_densities_sol = cell_densities.(sol.u)
        node_densities_sol = node_densities.(sol.u)
        for i in eachindex(sol)
            r = sol.u[i]
            q = cell_densities_sol[i]
            nq = node_densities_sol[i]
            for j in eachindex(r)
                if j == firstindex(r)
                    @test nq[j] ≈ 1 / (r[j+1] - r[j])
                    @test q[j] ≈ 1 / (r[j+1] - r[j])
                elseif j == lastindex(r)
                    @test nq[j] ≈ 1 / (r[j] - r[j-1])
                else
                    @test nq[j] ≈ 2 / (r[j+1] - r[j-1])
                    @test q[j] ≈ 1 / (r[j+1] - r[j])
                end
            end
        end

        # Compare with the continuum limit 
        let x = mb_prob.geometry.mesh_points
            fig = Figure(fontsize=36)
            colors = (:red, :blue, :darkgreen, :black, :orange, :magenta, :cyan, :yellow, :brown, :gray, :lightblue)
            ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"q(x, t)",
                title=L"(a): $q(x, t)$", titlealign=:left,
                width=600, height=300)
            [lines!(ax, sol.u[i], node_densities_sol[i], color=colors[i], linewidth=2) for i in eachindex(sol)]
            @views [lines!(ax, x .* mb_sol.u[i][end], mb_sol.u[i][begin:(end-1)], color=colors[i], linewidth=4, linestyle=:dashdot) for i in eachindex(sol)]
            resize_to_layout!(fig)
            fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
            @test_reference joinpath(fig_path, "step_function_moving_boundary.png") fig
            fig
        end

        # The steady problem should have diff = 0.2 
        diffs = diff(sol.u[end])
        @test all(x -> isapprox(x, 0.2, rtol=1e-1), diffs)
        sprob = SteadyCellProblem(prob)
        ssol = solve(sprob, DynamicSS(TRBDF2()))
        @test all(x -> isapprox(x, 0.2, rtol=1e-4), diff(ssol.u))
        @test ssol.u[end] ≈ 9.2 rtol = 1e-2
        @test sol.u[end][end] ≈ 9.2 rtol = 1e-1
    end
end

@testset "Proliferation with a Fixed Boundary" begin
    println("Starting Proliferation with a Fixed Boundary example.")
    force_law = (δ, p) -> p.k * (p.s - δ)
    force_law_parameters = (k=10.0, s=0.2)
    proliferation_law = (δ, p) -> max(zero(δ), p.β * p.K * (one(δ) - inv(p.K * δ)))
    proliferation_law_parameters = (β=1e-2, K=15.0)
    proliferation_period = 1e-2
    final_time = 50.0
    damping_constant = 1.0
    initial_condition = [LinRange(0, 15, 16); LinRange(15, 30, 32)] |> unique!
    prob = CellProblem(;
        force_law,
        force_law_parameters,
        proliferation_law,
        proliferation_law_parameters,
        proliferation_period,
        final_time,
        damping_constant,
        initial_condition)
    ens_prob = EnsembleProblem(prob)
    for sol in (solve(ens_prob, Tsit5(); trajectories=50, saveat=0.01), solve(ens_prob, Tsit5(); trajectories=50, saveat=0.01))
        sol = solve(ens_prob, Tsit5(); trajectories=50, saveat=0.01)

        fvm_prob = continuum_limit(prob, 1000; proliferation=true)
        fvm_sol = solve(fvm_prob, TRBDF2(linsolve=KLUFactorization()), saveat=0.01)

        q, r, means, lowers, uppers, knots = node_densities(sol)
        @inferred node_densities(sol)
        N, N_means, N_lowers, N_uppers = cell_numbers(sol)
        @inferred cell_numbers(sol)
        pde_N = map(eachindex(fvm_sol)) do i
            integrate_pde(fvm_sol, i)
        end
        @test all(≈(LinRange(0, 30, 500)), knots)

        # Test the statistics 
        for k in rand(1:length(sol), 20)
            for j in rand(1:length(sol[k]), 40)
                for i in rand(1:length(sol[k][j]), 60)
                    if i == 1
                        @test q[k][j][1] ≈ 1 / (r[k][j][2] - r[k][j][1])
                    elseif i == length(sol[k][j])
                        n = length(sol[k][j])
                        @test q[k][j][n] ≈ 1 / (r[k][j][n] - r[k][j][n-1])
                    else
                        @test q[k][j][i] ≈ 2 / (r[k][j][i+1] - r[k][j][i-1])
                    end
                    @test r[k][j][i] == sol[k][j][i]
                end
                @test N[k][j] ≈ length(sol[k][j]) - 1
            end
        end
        for j in rand(1:length(fvm_sol), 50)
            Nj = [N[k][j] for k in eachindex(sol)]
            @test mean(Nj) ≈ N_means[j]
            @test quantile(Nj, 0.025) ≈ N_lowers[j]
            @test quantile(Nj, 0.975) ≈ N_uppers[j]
            @test pde_N[j] ≈ DataInterpolations.integral(
                LinearInterpolation(fvm_sol.u[j], fvm_sol.prob.p.geometry.mesh_points),
                0.0, 30.0
            )
            for i in rand(1:length(knots[j]), 50)
                all_q = [LinearInterpolation(q[k][j], r[k][j])(knots[j][i]) for k in eachindex(sol)]
                @test mean(all_q) ≈ means[j][i]
                @test quantile(all_q, 0.025) ≈ lowers[j][i]
                @test quantile(all_q, 0.975) ≈ uppers[j][i]
            end
        end

        @test N_means[end] ≈ 450 rtol = 1e-2 # 30K

        fig = Figure(fontsize=33)
        colors = (:red, :blue, :darkgreen, :black, :magenta, :brown)
        plot_idx = (1, 1001, 2001, 3001, 4001, 5001)

        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"q(x, t)",
            title=L"(a): $q(x, t)$", titlealign=:left,
            width=600, height=300)
        [band!(ax, knots[i], lowers[i], uppers[i], color=(colors[j], 0.3)) for (j, i) in enumerate(plot_idx)]
        [lines!(ax, knots[i], means[i], color=colors[j], linewidth=2) for (j, i) in enumerate(plot_idx)]
        [lines!(ax, fvm_prob.geometry.mesh_points, fvm_sol.u[i], color=colors[j], linewidth=4, linestyle=:dashdot) for (j, i) in enumerate(plot_idx)]

        ax = Axis(fig[1, 2], xlabel=L"t", ylabel=L"N(t)",
            title=L"(b): $N(t)$", titlealign=:left,
            width=600, height=300)
        band!(ax, fvm_sol.t, N_lowers, N_uppers, color=(:blue, 0.3))
        lines!(ax, fvm_sol.t, N_means, color=:blue, linewidth=2)
        lines!(ax, fvm_sol.t, pde_N, color=:black, linewidth=4, linestyle=:dashdot)

        resize_to_layout!(fig)
        fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
        @test_reference joinpath(fig_path, "step_function_proliferation.png") fig by = psnr_equality(16.5)

        # Test the statistics when restricting to a specific set of simulation indices
        _indices = rand(eachindex(sol), 20)
        q, r, means, lowers, uppers, knots = node_densities(sol; indices=_indices)
        @inferred node_densities(sol; indices=_indices)
        N, N_means, N_lowers, N_uppers = cell_numbers(sol; indices=_indices)
        @inferred cell_numbers(sol; indices=_indices)
        @test all(≈(LinRange(0, 30, 500)), knots)
        for (enum_k, k) in enumerate(_indices)
            for j in rand(1:length(sol[k]), 40)
                for i in rand(1:length(sol[k][j]), 60)
                    if i == 1
                        @test q[enum_k][j][1] ≈ 1 / (r[enum_k][j][2] - r[enum_k][j][1])
                    elseif i == length(sol[k][j])
                        n = length(sol[k][j])
                        @test q[enum_k][j][n] ≈ 1 / (r[enum_k][j][n] - r[enum_k][j][n-1])
                    else
                        @test q[enum_k][j][i] ≈ 2 / (r[enum_k][j][i+1] - r[enum_k][j][i-1])
                    end
                    @test r[enum_k][j][i] == sol[k][j][i]
                end
                @test N[enum_k][j] ≈ length(sol[k][j]) - 1
            end
        end
        for j in rand(1:length(fvm_sol), 50)
            Nj = [N[k][j] for k in eachindex(_indices)]
            @test mean(Nj) ≈ N_means[j]
            @test quantile(Nj, 0.025) ≈ N_lowers[j]
            @test quantile(Nj, 0.975) ≈ N_uppers[j]
            @test pde_N[j] ≈ DataInterpolations.integral(
                LinearInterpolation(fvm_sol.u[j], fvm_sol.prob.p.geometry.mesh_points),
                0.0, 30.0
            )
            for i in rand(1:length(knots[j]), 50)
                all_q = [LinearInterpolation(q[k][j], r[k][j])(knots[j][i]) for k in eachindex(_indices)]
                @test mean(all_q) ≈ means[j][i]
                @test quantile(all_q, 0.025) ≈ lowers[j][i]
                @test quantile(all_q, 0.975) ≈ uppers[j][i]
            end
        end

        # Test the statistics with a specific interpolation function 
        _indices = eachindex(sol)
        q, r, means, lowers, uppers, knots = node_densities(sol; indices=_indices, interp_fnc=CubicSpline)
        @inferred node_densities(sol; indices=_indices, interp_fnc=CubicSpline)
        @test all(≈(LinRange(0, 30, 500)), knots)
        for (enum_k, k) in enumerate(_indices)
            for j in rand(1:length(sol[k]), 40)
                for i in rand(1:length(sol[k][j]), 60)
                    if i == 1
                        @test q[enum_k][j][1] ≈ 1 / (r[enum_k][j][2] - r[enum_k][j][1])
                    elseif i == length(sol[k][j])
                        n = length(sol[k][j])
                        @test q[enum_k][j][n] ≈ 1 / (r[enum_k][j][n] - r[enum_k][j][n-1])
                    else
                        @test q[enum_k][j][i] ≈ 2 / (r[enum_k][j][i+1] - r[enum_k][j][i-1])
                    end
                    @test r[enum_k][j][i] == sol[k][j][i]
                end
            end
        end
        for j in rand(1:length(fvm_sol), 50)
            for i in rand(1:length(knots[j]), 50)
                all_q = [CubicSpline(q[k][j], r[k][j])(knots[j][i]) for k in eachindex(_indices)]
                @test mean(all_q) ≈ means[j][i]
                @test quantile(all_q, 0.025) ≈ lowers[j][i]
                @test quantile(all_q, 0.975) ≈ uppers[j][i]
            end
        end

        # Using average leading edge 
        _indices = rand(eachindex(sol), 40)
        q, r, means, lowers, uppers, knots = node_densities(sol; indices=_indices, use_extrema=false)
        @inferred node_densities(sol; indices=_indices, use_extrema=false)
        @test all(≈(LinRange(0, 30, 500)), knots)
        for (enum_k, k) in enumerate(_indices)
            for j in rand(1:length(sol[k]), 40)
                for i in rand(1:length(sol[k][j]), 60)
                    if i == 1
                        @test q[enum_k][j][1] ≈ 1 / (r[enum_k][j][2] - r[enum_k][j][1])
                    elseif i == length(sol[k][j])
                        n = length(sol[k][j])
                        @test q[enum_k][j][n] ≈ 1 / (r[enum_k][j][n] - r[enum_k][j][n-1])
                    else
                        @test q[enum_k][j][i] ≈ 2 / (r[enum_k][j][i+1] - r[enum_k][j][i-1])
                    end
                    @test r[enum_k][j][i] == sol[k][j][i]
                end
            end
        end
        for j in rand(1:length(fvm_sol), 50)
            for i in rand(1:length(knots[j]), 50)
                all_q = [LinearInterpolation(q[k][j], r[k][j])(knots[j][i]) for k in eachindex(_indices)]
                @test mean(all_q) ≈ means[j][i]
                @test quantile(all_q, 0.025) ≈ lowers[j][i]
                @test quantile(all_q, 0.975) ≈ uppers[j][i]
            end
        end
    end
end

@testset "Proliferation with a Moving Boundary" begin
    println("Starting Proliferation with a Moving Boundary example.")
    force_law = (δ, p) -> p.k * (p.s - δ)
    force_law_parameters = (k=10.0, s=1)
    proliferation_law = (δ, p) -> max(zero(δ), p.β * p.K * (one(δ) - inv(p.K * δ)))
    proliferation_law_parameters = (β=1e-2, K=15.0)
    proliferation_period = 1e-2
    final_time = 50.0
    damping_constant = 1.0
    initial_condition = [LinRange(0, 15, 16); LinRange(15, 30, 32)] |> unique!
    prob = CellProblem(;
        force_law,
        force_law_parameters,
        proliferation_law,
        proliferation_law_parameters,
        proliferation_period,
        final_time,
        damping_constant,
        initial_condition,
        fix_right=false)
    for ens_prob in (EnsembleProblem(prob), EnsembleProblem(prob, sort=true))
        sol = solve(ens_prob, Tsit5(); trajectories=50, saveat=0.01)

        mb_prob = continuum_limit(prob, 1000; proliferation=true)
        mb_sol = solve(mb_prob, TRBDF2(linsolve=KLUFactorization()), saveat=0.01)

        q, r, means, lowers, uppers, knots = node_densities(sol)
        @inferred node_densities(sol)
        N, N_means, N_lowers, N_uppers = cell_numbers(sol)
        @inferred cell_numbers(sol)
        L, L_means, L_lowers, L_uppers = leading_edges(sol)
        pde_N = map(eachindex(mb_sol)) do i
            integrate_pde(mb_sol, i)
        end
        pde_L = map(mb_sol) do u
            u[end]
        end

        # Test the knots
        for j in eachindex(knots)
            a = Inf
            b = -Inf
            m = minimum(sol[k][j][begin] for k in eachindex(sol))
            M = maximum(sol[k][j][end] for k in eachindex(sol))
            @test knots[j] == LinRange(m, M, 500)
        end

        # Test the statistics
        for k in rand(1:length(sol), 20)
            for j in rand(1:length(sol[k]), 40)
                for i in rand(1:length(sol[k][j]), 60)
                    if i == 1
                        @test q[k][j][1] ≈ 1 / (r[k][j][2] - r[k][j][1])
                    elseif i == length(sol[k][j])
                        n = length(sol[k][j])
                        @test q[k][j][n] ≈ 1 / (r[k][j][n] - r[k][j][n-1])
                    else
                        @test q[k][j][i] ≈ 2 / (r[k][j][i+1] - r[k][j][i-1])
                    end
                    @test r[k][j][i] == sol[k][j][i]
                end
                @test N[k][j] ≈ length(sol[k][j]) - 1
            end
        end
        for j in rand(eachindex(mb_sol), 40)
            Nj = [N[k][j] for k in eachindex(sol)]
            @test @views mean(Nj) ≈ N_means[j]
            @test @views quantile(Nj, 0.025) ≈ N_lowers[j]
            @test @views quantile(Nj, 0.975) ≈ N_uppers[j]
            Lj = [L[k][j] for k in eachindex(sol)]
            @test @views mean(Lj) ≈ L_means[j]
            @test @views quantile(Lj, 0.025) ≈ L_lowers[j]
            @test @views quantile(Lj, 0.975) ≈ L_uppers[j]
            @test pde_N[j] ≈ DataInterpolations.integral(
                LinearInterpolation(mb_sol.u[j][begin:(end-1)], mb_sol.u[j][end] * mb_sol.prob.p.geometry.mesh_points),
                0.0, mb_sol.u[j][end]
            )
            @test pde_L[j] ≈ mb_sol.u[j][end]
            for i in rand(eachindex(knots[j]), 60)
                all_q = max.(0.0, [LinearInterpolation(q[k][j], r[k][j])(knots[j][i]) * (knots[j][i] ≤ r[k][j][end]) for k in eachindex(sol)])
                @test mean(all_q) ≈ means[j][i] rtol = 1e-3
                @test quantile(all_q, 0.025) ≈ lowers[j][i] rtol = 1e-3
                @test quantile(all_q, 0.975) ≈ uppers[j][i] rtol = 1e-3
            end
        end

        fig = Figure(fontsize=33)
        colors = (:red, :blue, :darkgreen, :black, :magenta, :brown)
        plot_idx = (1, 1001, 2001, 3001, 4001, 5001)

        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"q(x, t)",
            title=L"(a): $q(x, t)$", titlealign=:left,
            width=600, height=300)
        [band!(ax, knots[i], lowers[i], uppers[i], color=(colors[j], 0.3)) for (j, i) in enumerate(plot_idx)]
        [lines!(ax, knots[i], means[i], color=colors[j], linewidth=2) for (j, i) in enumerate(plot_idx)]
        @views [lines!(ax, mb_sol.u[i][end] .* mb_prob.geometry.mesh_points, mb_sol.u[i][begin:(end-1)], color=colors[j], linewidth=4, linestyle=:dashdot) for (j, i) in enumerate(plot_idx)]

        ax = Axis(fig[1, 2], xlabel=L"t", ylabel=L"N(t)",
            title=L"(b): $N(t)$", titlealign=:left,
            width=600, height=300)
        band!(ax, mb_sol.t, N_lowers, N_uppers, color=(:blue, 0.3))
        lines!(ax, mb_sol.t, N_means, color=:blue, linewidth=2)
        lines!(ax, mb_sol.t, pde_N, color=:black, linewidth=4, linestyle=:dashdot)

        ax = Axis(fig[1, 3], xlabel=L"t", ylabel=L"L(t)",
            title=L"(c): $L(t)$", titlealign=:left,
            width=600, height=300)
        band!(ax, mb_sol.t, L_lowers, L_uppers, color=(:blue, 0.3))
        lines!(ax, mb_sol.t, L_means, color=:blue, linewidth=2)
        lines!(ax, mb_sol.t, pde_L, color=:black, linewidth=4, linestyle=:dashdot)

        resize_to_layout!(fig)
        fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
        @test_reference joinpath(fig_path, "step_function_proliferation_moving_boundary.png") fig by = psnr_equality(15)

        # Test the statistics when restricting to a specific set of simulation indices
        _indices = rand(eachindex(sol), 20)
        q, r, means, lowers, uppers, knots = node_densities(sol; indices=_indices)
        @inferred node_densities(sol; indices=_indices)
        N, N_means, N_lowers, N_uppers = cell_numbers(sol; indices=_indices)
        @inferred cell_numbers(sol; indices=_indices)
        L, L_means, L_lowers, L_uppers = leading_edges(sol; indices=_indices)
        for j in eachindex(knots)
            a = Inf
            b = -Inf
            m = minimum(sol[k][j][begin] for k in _indices)
            M = maximum(sol[k][j][end] for k in _indices)
            @test knots[j] == LinRange(m, M, 500)
        end
        for (enum_k, k) in enumerate(_indices)
            for j in rand(1:length(sol[k]), 40)
                for i in rand(1:length(sol[k][j]), 60)
                    if i == 1
                        @test q[enum_k][j][1] ≈ 1 / (r[enum_k][j][2] - r[enum_k][j][1])
                    elseif i == length(sol[k][j])
                        n = length(sol[k][j])
                        @test q[enum_k][j][n] ≈ 1 / (r[enum_k][j][n] - r[enum_k][j][n-1])
                    else
                        @test q[enum_k][j][i] ≈ 2 / (r[enum_k][j][i+1] - r[enum_k][j][i-1])
                    end
                    @test r[enum_k][j][i] == sol[k][j][i]
                end
                @test N[enum_k][j] ≈ length(sol[k][j]) - 1
            end
        end
        for j in rand(eachindex(mb_sol), 40)
            Nj = [N[k][j] for k in eachindex(_indices)]
            @test @views mean(Nj) ≈ N_means[j]
            @test @views quantile(Nj, 0.025) ≈ N_lowers[j]
            @test @views quantile(Nj, 0.975) ≈ N_uppers[j]
            Lj = [L[k][j] for k in eachindex(_indices)]
            @test @views mean(Lj) ≈ L_means[j]
            @test @views quantile(Lj, 0.025) ≈ L_lowers[j]
            @test @views quantile(Lj, 0.975) ≈ L_uppers[j]
            @test pde_N[j] ≈ DataInterpolations.integral(
                LinearInterpolation(mb_sol.u[j][begin:(end-1)], mb_sol.u[j][end] * mb_sol.prob.p.geometry.mesh_points),
                0.0, mb_sol.u[j][end]
            )
            @test pde_L[j] ≈ mb_sol.u[j][end]
            for i in rand(eachindex(knots[j]), 60)
                all_q = max.(0.0, [LinearInterpolation(q[k][j], r[k][j])(knots[j][i]) * (knots[j][i] ≤ r[k][j][end]) for k in eachindex(_indices)])
                @test mean(all_q) ≈ means[j][i] rtol = 1e-3
                @test quantile(all_q, 0.025) ≈ lowers[j][i] rtol = 1e-3
                @test quantile(all_q, 0.975) ≈ uppers[j][i] rtol = 1e-3
            end
        end

        # Test the statistics with a specific interpolation function 
        _indices = rand(eachindex(sol), 20)
        q, r, means, lowers, uppers, knots = node_densities(sol; indices=_indices, interp_fnc=CubicSpline)
        @inferred node_densities(sol; indices=_indices, interp_fnc=CubicSpline)
        for j in eachindex(knots)
            a = Inf
            b = -Inf
            m = minimum(sol[k][j][begin] for k in _indices)
            M = maximum(sol[k][j][end] for k in _indices)
            @test knots[j] == LinRange(m, M, 500)
        end
        for (enum_k, k) in enumerate(_indices)
            for j in rand(1:length(sol[k]), 40)
                for i in rand(1:length(sol[k][j]), 60)
                    if i == 1
                        @test q[enum_k][j][1] ≈ 1 / (r[enum_k][j][2] - r[enum_k][j][1])
                    elseif i == length(sol[k][j])
                        n = length(sol[k][j])
                        @test q[enum_k][j][n] ≈ 1 / (r[enum_k][j][n] - r[enum_k][j][n-1])
                    else
                        @test q[enum_k][j][i] ≈ 2 / (r[enum_k][j][i+1] - r[enum_k][j][i-1])
                    end
                    @test r[enum_k][j][i] == sol[k][j][i]
                end
            end
        end
        for j in rand(eachindex(mb_sol), 40)
            for i in rand(eachindex(knots[j]), 60)
                all_q = max.(0.0, [CubicSpline(q[k][j], r[k][j])(knots[j][i]) * (knots[j][i] ≤ r[k][j][end]) for k in eachindex(_indices)])
                @test mean(all_q) ≈ means[j][i] rtol = 1e-3
                @test quantile(all_q, 0.025) ≈ lowers[j][i] rtol = 1e-3
                @test quantile(all_q, 0.975) ≈ uppers[j][i] rtol = 1e-3
            end
        end

        # Using the average leading edge
        L = leading_edges(sol).L
        _L = stack(L)
        _indices = rand(eachindex(sol), 20)
        _L = _L[:, _indices]
        _mL = mean.(eachrow(_L))
        q, r, means, lowers, uppers, knots = node_densities(sol; indices=_indices, use_extrema=false)
        @inferred node_densities(sol; indices=_indices, use_extrema=false)
        for j in eachindex(knots)
            a = mean(sol[k][j][begin] for k in _indices)
            b = mean(sol[k][j][end] for k in _indices)
            @test knots[j] ≈ LinRange(a, b, 500)
            @test knots[j][end] ≈ _mL[j]
        end
        for (enum_k, k) in enumerate(_indices)
            for j in rand(1:length(sol[k]), 40)
                for i in 1:length(sol[k][j])
                    if i == 1
                        @test q[enum_k][j][1] ≈ 1 / (r[enum_k][j][2] - r[enum_k][j][1])
                    elseif i == length(sol[k][j])
                        n = length(sol[k][j])
                        @test q[enum_k][j][n] ≈ 1 / (r[enum_k][j][n] - r[enum_k][j][n-1])
                    else
                        @test q[enum_k][j][i] ≈ 2 / (r[enum_k][j][i+1] - r[enum_k][j][i-1])
                    end
                    @test r[enum_k][j][i] == sol[k][j][i]
                end
            end
        end
        for j in rand(eachindex(mb_sol), 40)
            for i in eachindex(knots[j])
                all_q = max.(0.0, [LinearInterpolation(q[k][j], r[k][j])(knots[j][i]) * (knots[j][i] ≤ r[k][j][end]) for k in eachindex(_indices)])
                @test mean(all_q) ≈ means[j][i] rtol = 1e-3
                @test quantile(all_q, 0.025) ≈ lowers[j][i] rtol = 1e-3
                @test quantile(all_q, 0.975) ≈ uppers[j][i] rtol = 1e-3
            end
        end
    end
end