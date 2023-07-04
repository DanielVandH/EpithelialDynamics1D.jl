var documenterSearchIndex = {"docs":
[{"location":"examples/","page":"Examples","title":"Examples","text":"CurrentModule = EpithelialDynamics1D","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This section gives some examples for how the package can be used. We consider the same type of initial condition, migration, and prloiferation mechanisms, but break the example into four parts:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Fixed boundaries, no proliferation;\nFixed boundaries, proliferation;\nFree boundaries, no proliferation;\nFree boundaries, proliferation.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The initial condition we use is a step function density on 0 leq x leq 30, with cells in 0 15 of low density and cells in 15 30 of great density:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"initial_condition = [LinRange(0, 15, 16); LinRange(15, 30, 32)] |> unique!\nfig = Figure(fontsize=33)\nax = Axis(fig[1, 1], xlabel=L\"x\",width=600,height=200)\nscatter!(ax, initial_condition, zero(initial_condition),color=:black,markersize=13)\nhideydecorations!(ax)\nresize_to_layout!(fig)\nfig","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<figure>\n    <img src='../figures/step_function_initial_condition.png', alt'Step function initial condition'><br>\n</figure>","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"When considering cell migration, the force law we use is the linear law F(delta) = k(s - delta), where k = 10 is the spring constant in each example, and s = 02 is the resting spring length except in the last example where s=1. When including proliferation, we use the logistic law G(delta) = max0 beta K1 - 1(Kdelta), where beta = 001 is the intrinsic proliferation rate and K=15 is the cell carrying capacity density. This choice of proliferation law slows down the growth of a cell population when there are many cells packed together, ensuring that a steady state can be reached (when there is no moving boundary). (To see that this is a logistic law, note that in the continuum limit the corresponding reaction term is qG(1q) = max0 beta K q1 - qK, which is basically the same term in the Fisher equation.)","category":"page"},{"location":"examples/#Example-I:-Fixed-boundaries,-no-proliferation","page":"Examples","title":"Example I: Fixed boundaries, no proliferation","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"We start with a problem that fixes both boundaries and only includes cell migration. The first step we take is to define the problem itself:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using EpithelialDynamics1D \nforce_law = (δ, p) -> p.k * (p.s - δ)\nforce_law_parameters = (k=10.0, s=0.2)\nfinal_time = 100.0\ndamping_constant = 1.0\ninitial_condition = [LinRange(0, 15, 16); LinRange(15, 30,32)] |> unique!\nprob = CellProblem(;\n    force_law,\n    force_law_parameters,\n    final_time,\n    damping_constant,\n    initial_condition)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This problem can be solved like any other problem from DifferentialEquations.jl. I find that Tsit5() is typically the fastest:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using OrdinaryDiffEq\nsol = solve(prob, Tsit5(), saveat=10.0)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Let's compare this solution to its continuum limit. From Baker et al. (2019), the continuum limit is given by ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign*\nbeginarrayrcll\ndfracpartial qpartial t  =  dfracpartialpartial xleft(D(q)dfracpartial qpartial xright)  0  x  Lt09pt\ndfracpartial qpartial x  =  0  x = 0t09pt\ndfracpartial qpartial x  =  0  x = Lt09pt\nq(x 0)  =  q_0(x)  0 leq x leq L\nendarray\nendalign*","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"These PDEs are solved using FiniteVolumeMethod1D.jl. This dependent variable q is the cell density.  There are two possible ways to define densities:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"cell_densities: The cell_densities function defines the density of a cell (x_i x_i+1) as the reciprocal length q_i = 1(x_i+1 - x_i).\nnode_densities: The node_densities function, the one that actually gets used in the continuum limit, assigns densities to nodes x_i rather than cells (x_i x_i+1). The density at a node x_i is defined by q_i = 2(x_i+1 - x_i-1) if 1  i  n, where n is the number of cells at the given time, or q_1 = 1(x_2 - x_1) and q_n = 1(x_n - x_n-1) at the endpoints.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The latter definition is used in the continuum limit - the initial condition q_0(x) is a piecewise linear interpolant through the cell densities at the initial time from the cell problem's initial condition. The diffusion function D(q) is given by ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"D(q) = -dfrac1eta q^2Fleft(dfrac1qright)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"which for our force law is D(q) = alphaq^2, alpha = keta.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The function continuum_limit constructs the FVMProblem defining this continuum limit:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"fvm_prob = continuum_limit(prob, 1000)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The second argument 1000 defines the number of mesh points to use for the continuum limit, using the piecewise linear interpolant q_0(x) to compute the densities at each mesh point. Now let's solve it:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using LinearSolve\nfvm_sol = solve(fvm_prob, TRBDF2(linsolve=KLUFactorization()), saveat=10.0)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"If the simulation is correct, noting that our value of k is sufficiently high so that the continuum limit should actually work, then fvm_sol should be a good match to the densities from sol. To verify this, let us make a plot:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"node_densities_sol = node_densities.(sol.u)\nusing CairoMakie\nlet x = fvm_prob.geometry.mesh_points\n    fig = Figure(fontsize=36)\n    colors = (:red, :blue, :darkgreen, :black, :orange, :magenta,:cyan, :yellow, :brown, :gray, :lightblue)\n    ax = Axis(fig[1, 1], xlabel=L\"x\", ylabel=L\"q(x, t)\",\n        title=L\"(a): $q(x, t)$\", titlealign=:left,\n        width=600, height=300)\n    [lines!(ax, sol.u[i], node_densities_sol[i], color=colors[i],linewidth=2) for i in eachindex(sol)]\n    [lines!(ax, x, fvm_sol.u[i], color=colors[i], linewidth=4,linestyle=:dashdot) for i in eachindex(sol)]\n    resize_to_layout!(fig)\n    fig\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<figure>\n    <img src='../figures/step_function.png', alt'Solution compared to continuum limit'><br>\n</figure>","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We see that the continuum limit lines up perfectly with the discrete results, and the solution approaches the steady state that is the average of the two densities from the initial condition. We provide a SteadyCellProblem for computing this steady state without running any simulation:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using SteadyStateDiffEq\nsprob = SteadyCellProblem(prob)\nssol = solve(sprob, DynamicSS(TRBDF2()))\nq = cell_densities(ssol.u)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> q = cell_densities(ssol.u)\n46-element Vector{Float64}:\n 1.5333521372657852\n 1.533352049592097\n 1.5333518746535277\n 1.5333516132657792\n 1.5333512666476505\n ⋮\n 1.533315400438977\n 1.5333150538372589\n 1.5333147924620867\n 1.5333146175320538\n 1.5333145298626596","category":"page"},{"location":"examples/#Example-II:-Moving-boundary,-no-proliferation","page":"Examples","title":"Example II: Moving boundary, no proliferation","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now let's allow the rightmost boundary to be free. This is done by setting fix_right = false in the CellProblem constructor.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using EpithelialDynamics1D, OrdinaryDiffEq\n\nforce_law = (δ, p) -> p.k * (p.s - δ)\nforce_law_parameters = (k=10.0, s=0.2)\nfinal_time = 500.0\ndamping_constant = 1.0\ninitial_condition = [LinRange(0, 15, 16); LinRange(15, 30, 32)] |>unique!\nprob = CellProblem(;\n    force_law,\n    force_law_parameters,\n    final_time,\n    damping_constant,\n    initial_condition,\n    fix_right=false)\n\nsol = solve(prob, Tsit5(), saveat=50.0)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The continuum limit that we can compare to in this case is this given by, letting L(t) denote the position of the rightmost node (the leading edge) at the time t.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign*\nbeginarrayrcll\ndfracpartial qpartial t  =  dfracpartialpartial xleft(D(q)dfracpartial qpartial xright)  0  x  L(t)t09pt\ndfracpartial qpartial x  =  0  x = 0 9pt\ndfrac1etaFleft(dfrac1qright) + dfracD(q)2qdfracpartial qpartial x  =  0  x=L(t)t09pt\ndfracmathrm dLmathrm dt  =  -dfrac1qD(q)dfracpartial qpartial x  x = L(t)t0 9pt\nq(x 0)  =  q_0(x)  0 leq x leq L(0) 9pt\nL(0)  =  L_0\nendarray\nendalign*","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"These PDEs are solved using MovingBoundaryProblems1D.jl. The initial endpoint L_0 is given by the position of the rightmost node at t=0, in this case L_0 = 30. Once again, we construct the corresponding problem, in this case an MBProblem, using continuum_limit:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using LinearSolve\n\nmb_prob = continuum_limit(prob, 1000)\nmb_sol = solve(mb_prob, TRBDF2(linsolve=KLUFactorization()), saveat=50.0)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The following code then plots and compares the two solutions.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"node_densities_sol = node_densities.(sol.u)\nusing CairoMakie\nlet x = mb_prob.geometry.mesh_points\n    fig = Figure(fontsize=36)\n    colors = (:red, :blue, :darkgreen, :black, :orange, :magenta,:cyan, :yellow, :brown, :gray, :lightblue)\n    ax = Axis(fig[1, 1], xlabel=L\"x\", ylabel=L\"q(x, t)\",\n        title=L\"(a): $q(x, t)$\", titlealign=:left,\n        width=600, height=300)\n    [lines!(ax, sol.u[i], node_densities_sol[i], color=colors[i],linewidth=2) for i in eachindex(sol)]\n    @views [lines!(ax, x .* mb_sol.u[i][end], mb_sol.u[i][begin(end-1)], color=colors[i], linewidth=4, linestyle=:dashdot)for i in eachindex(sol)]\n    resize_to_layout!(fig)\n    fig\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<figure>\n    <img src='../figures/step_function_moving_boundary.png', alt'Solution compared to continuum limit'><br>\n</figure>","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We see that the solutions match. The cells appear to be retreating from L(0) = 30, contracting inwards until a steady state is eventually reached. In this steady state,  the cells seem to pack into the interval 0 10. This makes sense - we start with 46 cells:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> length(initial_condition) - 1\n46","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"With a resting spring length of s = 02, we can fit 46 cells into the interval 0 46 cdot 02 = 0 92. We can verify the properties of this steady state by either looking at the last time from the simulation, or by computing the steady state directly:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> diffs = diff(sol.u[end])\n46-element Vector{Float64}:\n 0.20233969090471757\n 0.20203667605421724\n 0.20263100257524164\n 0.20172744674985177\n 0.20290519447397481\n ⋮\n 0.19964286225916084\n 0.20089946315869867\n 0.1997809643420574\n 0.20045281417783656\n 0.19992619110920273","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Indeed, the length of each cell is approximately s = 02, and the endpoint at the last time is:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> sol.u[end][end]\n9.267001076147542","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"which is approaching 92. If we look at the steady state:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using SteadyStateDiffEq \nsprob = SteadyCellProblem(prob)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> ssol = solve(sprob, DynamicSS(TRBDF2()))\nu: 47-element Vector{Float64}:\n -3.8321946143843015e-20\n  0.20000125473810257\n  0.4000025080445262\n  0.6000037584892256\n  0.8000050046454205\n  ⋮\n  8.400036722439138\n  8.600036891490856\n  8.800037018448712\n  9.000037103167848\n  9.200037145551596","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We see that the cells all have lengths approximately equal to s, and  lim_t to infty L(t) = 92 as predicted.","category":"page"},{"location":"examples/#Example-III:-Fixed-boundary,-with-proliferation","page":"Examples","title":"Example III: Fixed boundary, with proliferation","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now let us go back to the fixed boundary problem, but now include proliferation. Remember that the proliferation law we use is the logistic law G(delta) = max0 beta K1 - 1(Kdelta).","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The procedure for solving problems with proliferation is different than without. Since the proliferation mechanism is stochastic, we need to simulate the system many times to capture the average behaviour. We use the ensemble solution features from DifferentialEquations.jl to do this, using the trajectories keyword to specify how many simulations of the systems we want.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using EpithelialDynamics1D, OrdinaryDiffEq \n\nforce_law = (δ, p) -> p.k * (p.s - δ)\nforce_law_parameters = (k=10.0, s=0.2)\nproliferation_law = (δ, p) -> max(zero(δ), p.β * p.K * (one(δ) -inv(p.K * δ)))\nproliferation_law_parameters = (β=1e-2, K=15.0)\nproliferation_period = 1e-2\nfinal_time = 50.0\ndamping_constant = 1.0\ninitial_condition = [LinRange(0, 15, 16); LinRange(15, 30, 32)] |>unique!\nprob = CellProblem(;\n    force_law,\n    force_law_parameters,\n    proliferation_law,\n    proliferation_law_parameters,\n    proliferation_period,\n    final_time,\n    damping_constant,\n    initial_condition)\nens_prob = EnsembleProblem(prob)\nsol = solve(ens_prob, Tsit5(); trajectories=50, saveat=0.01)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The continuum limit for this problem is similar to the problem without proliferation, except the PDE has a reaction term:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"dfracpartial qpartial t = dfracpartialpartial xleft(D(q)dfracpartial qpartial xright) + qGleft(dfrac1qright)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"the other boundary and initial conditions are the same as in the first example. As in the first example, PDEs of this form are solved using FiniteVolumeMethod1D.jl.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using LinearSolve \n\nfvm_prob = continuum_limit(prob, 1000; proliferation=true)\nfvm_sol = solve(fvm_prob, TRBDF2(linsolve=KLUFactorization()), saveat=0.01)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Note that the proliferation=true keyword argument is needed in continuum_limit, else the continuum limit without proliferation is returned. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We provide several functions for computing statistics from the EnsembleSolution, sol. The function node_densities returns a NamedTuple with the following properties:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"q: This is a vector-of-vector-of-vectors, where q[k][j][i] is the density of the ith node from the jth timepoint of the kth simulation.\nr: This is a vector-of-vector-of-vectors, where r[k][j][i] is the position of the ith node from the jth timepoint of the kth simulation.\nknots: To summarise the behaviour of the system at each time, we define a grid of knots at each time (defaults to 500 knots for each time). These knots are used to evaluate the piecewise linear interpolant of the densities at each time for each simulation, allowing us to summarise the densities at common knots for each time. With this property, knots[j] is the set of knots used for the jth time, where knots[j][begin] is the minimum of all cell positions from each simulation at the jth time, and knots[j][end] is the corresponding maximum. In this case, the minimum and maximum for each time are just 0 and 30, respectively.\nmeans: This is a vector-of-vectors, where means[j] is a vector of average densities at each knot in knots[j].\nlowers: This is a vector-of-vectors, where lowers[j] are the lower limits of the confidence intervals for the densities at each knot in knots[j]. The significance level of this confidence interval is alpha=005 by default, meaning these are (100alpha2) = 25 quantiles.\nuppers: Similar to lowers, except these are the upper limits of the corresponding confidence intervals, i.e. the 100(1-alpha2) = 975 ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Another function that we provide is cell_numbers, used for obtaining cell numbers at each time and summarising for each simulation. The returned result is a NamedTuple with the following properties:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"N: This is a vector-of-vectors, with N[k][j] the number of cells at the jth time of the kth simulation.\nmeans: This is a vector, with N[j] the average number of cells at the jth time.\nlowers: This is a vector, with N[j] the lower limit of the confidence interval of the cell numbers at the jth time. The significance level of this confidence interval is alpha=005 by default, meaning these are (100alpha2) = 25 quantiles.\nuppers: Similar to lowers, except these are the upper limits of the corresponding confidence intervals, i.e. the 100(1-alpha2) = 975.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Another useful function for comparing with the PDEs is integrate_pde, which integrates the PDEs at a given time, returning N(t) = int_0^L(t) q(x t) mathrm dx, an estimate for the number of cells at the time t. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Let's now compute these statistics:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(; q, r, means, lowers, uppers,knots) = node_densities(sol)\nN, N_means, N_lowers, N_uppers =cell_numbers(sol)\npde_N = map(eachindex(fvm_sol)) do i\n    integrate_pde(fvm_sol, i)\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now we can plot.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using CairoMakie \n\nfig = Figure(fontsize=33)\ncolors = (:red, :blue, :darkgreen, :black, :magenta, :brown)\nplot_idx = (1, 1001, 2001, 3001, 4001, 5001)\nax = Axis(fig[1, 1], xlabel=L\"x\",ylabel=L\"q(x, t)\",\n    title=L\"(a): $q(x, t)$\",titlealign=:left,\n    width=600, height=300)\n[band!(ax, knots[i], lowers[i],uppers[i], color=(colors[j], 0.3))for (j, i) in enumerate(plot_idx)]\n[lines!(ax, knots[i], means[i],color=colors[j], linewidth=2) for (j, i) in enumerate(plot_idx)]\n[lines!(ax, fvm_prob.geometrymesh_points, fvm_sol.u[i],color=colors[j], linewidth=4,linestyle=:dashdot) for (j, i) in enumerate(plot_idx)]\nax = Axis(fig[1, 2], xlabel=L\"t\",ylabel=L\"N(t)\",\n    title=L\"(b): $N(t)$\",titlealign=:left,\n    width=600, height=300)\nband!(ax, fvm_sol.t, N_lowers,N_uppers, color=(:blue, 0.3))\nlines!(ax, fvm_sol.t, N_means,color=:blue, linewidth=2)\nlines!(ax, fvm_sol.t, pde_N,color=:black, linewidth=4,linestyle=:dashdot)\nresize_to_layout!(fig)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<figure>\n    <img src='../figures/step_function_proliferation.png', alt'Solution compared to continuum limit'><br>\n</figure>","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The confidence intervals are shown in the shaded regions, with the mean values given by a solid curve. We see that the continuum limit is a good match for the mean behaviour of both q(x t) and N(t). The density q(x t) reaches a steady state with q(x t) to 15 as t to infty for each x, which we expect given the cell carrying capacity density K = 15. If we had used e.g. a constant proliferation law G(delta) = beta, we would see growth indefinitely, so the logistic law is especially nice for this reason. The cell numbers reach a limit around 450, which makes sense since int_0^30 Kmathrm dx = 30K = 450.","category":"page"},{"location":"examples/#Example-IV:-Moving-boundary,-with-proliferation","page":"Examples","title":"Example IV: Moving boundary, with proliferation","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now let us consider proliferation with a moving boundary.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using EpithelialDynamics1D, OrdinaryDiffEq\nforce_law = (δ, p) -> p.k * (p.s - δ)\nforce_law_parameters = (k=10.0, s=1)\nproliferation_law = (δ, p) -> max(zero(δ), p.β * p.K *(one(δ) - inv(p.K * δ)))\nproliferation_law_parameters = (β=1e-2, K=15.0)\nproliferation_period = 1e-2\nfinal_time = 30.0\ndamping_constant = 1.0\ninitial_condition = [LinRange(0, 15, 16); LinRange(15,30, 32)] |> unique!\nprob = CellProblem(;\n    force_law,\n    force_law_parameters,\n    proliferation_law,\n    proliferation_law_parameters,\n    proliferation_period,\n    final_time,\n    damping_constant,\n    initial_condition,\n    fix_right=false)\nens_prob = EnsembleProblem(prob)\nsol = solve(ens_prob, Tsit5(); trajectories=50, saveat=001)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The continuum limit is the same as it was in Example II, except now the PDE is the one given in Example III. We construct this continuum limit as follows, noting again that MovingBoundaryProblems1D.jl solves the moving boundary problem:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using LinearSolve\nmb_prob = continuum_limit(prob, 1000; proliferation=true)\nmb_sol = solve(mb_prob, TRBDF2(linsolve=KLUFactorization()), saveat=0.01)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"As in Example III, we can compute statistics from our ensemble solutions, as well as the corresponding values from the continuum limit:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(; q, r, means, lowers, uppers, knots) =node_densities(sol)\n@inferred node_densities(sol)\nN, N_means, N_lowers, N_uppers = cell_numbers(sol)\n@inferred cell_numbers(sol)\nL, L_means, L_lowers, L_uppers = leading_edge(sol)\npde_N = map(eachindex(mb_sol)) do i\n    integrate_pde(mb_sol, i)\nend\npde_L = map(mb_sol) do u\n    u[end]\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Let's now plot these results.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using CairoMakie\nfig = Figure(fontsize=33)\ncolors = (:red, :blue, :darkgreen,:black, :magenta, :brown)\nplot_idx = (1, 1001, 2001, 3001,4001, 5001)\nax = Axis(fig[1, 1], xlabel=L\"x\",ylabel=L\"q(x, t)\",\n    title=L\"(a): $q(x, t)$\",titlealign=:left,\n    width=600, height=300)\n[band!(ax, knots[i], lowers[i],uppers[i], color=(colors[j], 0.3))for (j, i) in enumerate(plot_idx)]\n[lines!(ax, knots[i], means[i],color=colors[j], linewidth=2) for(j, i) in enumerate(plot_idx)]\n@views [lines!(ax, mb_sol.u[i][end] * mb_prob.geometry.mesh_points,mb_sol.u[i][begin:(end-1)],color=colors[j], linewidth=4,linestyle=:dashdot) for (j, i) in enumerate(plot_idx)]\nax = Axis(fig[1, 2], xlabel=L\"t\",ylabel=L\"N(t)\",\n    title=L\"(b): $N(t)$\",titlealign=:left,\n    width=600, height=300)\nband!(ax, mb_sol.t, N_lowers,N_uppers, color=(:blue, 0.3))\nlines!(ax, mb_sol.t, N_means,color=:blue, linewidth=2)\nlines!(ax, mb_sol.t, pde_N,color=:black, linewidth=4,linestyle=:dashdot)\nax = Axis(fig[1, 3], xlabel=L\"t\",ylabel=L\"L(t)\",\n    title=L\"(c): $L(t)$\",titlealign=:left,\n    width=600, height=300)\nband!(ax, mb_sol.t, L_lowers,L_uppers, color=(:blue, 0.3))\nlines!(ax, mb_sol.t, L_means,color=:blue, linewidth=2)\nlines!(ax, mb_sol.t, pde_L,color=:black, linewidth=4,linestyle=:dashdot)\nresize_to_layout!(fig)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<figure>\n    <img src='../figures/step_function_proliferation_moving_boundary.png', alt'Solution compared to continuum limit'><br>\n</figure>","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Once again, the continuum limit is a great match.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = EpithelialDynamics1D","category":"page"},{"location":"#EpithelialDynamics1D","page":"Home","title":"EpithelialDynamics1D","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for EpithelialDynamics1D.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is a package for simulating epithelial dynamics in one dimension, supporting cell migration and cell proliferation, implementing the model in Baker et al. (2019). In this model, cells are represented as intervals between points, with the ith cell given by (x_i x_i+1), i=1ldotsn-1. The node positions x_i are governed by the differential equation","category":"page"},{"location":"","page":"Home","title":"Home","text":"etadfracmathrm dx_imathrm dt = Fleft(x_i - x_i-1right) - Fleft(x_i+1 - x_iright) quad i=2ldotsn-1","category":"page"},{"location":"","page":"Home","title":"Home","text":"assuming x_1  x_2  cdots  x_n. If the left boundary is fixed, then mathrm dx_1mathrm dt = 0, otherwise etamathrm dx_1mathrm dt = -F(x_2 - x_1). Similarly, if the right boundary is fixed then mathrm dx_nmathrm dt = 0, otherwise etamathrm dx_nmathrm dt = F(x_n - x_n-1). The parameter eta is called the damping constant or the viscosity coefficient and F is the force law. In this model, cells are modelled as springs, whose forces are governed by this force law F. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The proliferation mechanism that we use assumes that at most one cell can divide at a time. Let C_i(t) denote the event that the ith cell divides in the time interval t t+mathrm dt). We assume that PrC_i(t) = G_imathrm dt, where G_i = G(x_i+1 - x_i) for some proliferation law G. When this division occurs we place a new position at (x_i + x_i+1)2, adjusting the indices of the x_i accordingly so that they remain sorted. We implement this mechanism by attempting a proliferation event every mathrm dt = Delta t units of time, so that the probability of a proliferation event occuring at the time t is sum_i=1^n-1 G_iDelta t, and the probability that the ith cell proliferates, given that a proliferation event does occur, is G_isum_j=1^n-1 G_j. DiffEqCallbacks.jl is used to implement the callback that attempts this event periodically while simulating the problem.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Examples of how to simulate these systems are given in the sidebar.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [EpithelialDynamics1D]","category":"page"},{"location":"#EpithelialDynamics1D.CellProblem","page":"Home","title":"EpithelialDynamics1D.CellProblem","text":"CellProblem{F,DF,DP,G,GP}\n\nStruct for representing a cell simulation problem. \n\nFields\n\nforce_law::F: The force law, e.g. (δ, p) -> p.k * (p.s - δ).\nforce_law_parameters::Fp: The parameters for the force law.\nproliferation_law::G = (δ, p) -> zero(δ): The proliferation law, e.g. (δ, p) -> p.β.\nproliferation_law_parameters::Gp = nothing: The parameters for the proliferation law.\nproliferation_period::Float64 = 1e-2: How often to attempt a proliferation event.\ninitial_time::Float64 = 0.0: The initial time.\nfinal_time::Float64: The final time.\nfix_left::Bool = true: Whether to fix the left endpoint.\nfix_right::Bool = true: Whether to fix the right endpoint.\ndamping_constant::Float64: The damping constant.\ninitial_condition::Vector{Float64}: The initial cell positions. Must be sorted.\n\nSolving\n\nYou can solve a CellProblem like you would solve a problem from DifferentialEquations.jl, e.g. with \n\nsolve(prob, Tsit5(), saveat = 0.1)\n\n\n\n\n\n","category":"type"},{"location":"#EpithelialDynamics1D.SteadyCellProblem","page":"Home","title":"EpithelialDynamics1D.SteadyCellProblem","text":"SteadyCellProblem{M}\n\nDefines a steady-state problem for a CellProblem. Only has a  single field, prob::CellProblem.\n\nYou can solve this problem as you would a NonlinearProblem, e.g. with \n\nsolve(prob, alg)\n\nwhere alg is e.g. DynamicSS(TRBDF2())) from SteadyStateDiffEq.jl.\n\n\n\n\n\n","category":"type"},{"location":"#FiniteVolumeMethod1D.FVMProblem","page":"Home","title":"FiniteVolumeMethod1D.FVMProblem","text":"FVMProblem(prob::CellProblem,\n    mesh_points=copy(prob.initial_condition);\n    diffusion_function=:continuum,\n    diffusion_parameters=nothing,\n    reaction_function=:continuum,\n    reaction_parameters=nothing,\n    proliferation=false)\n\nConstructs an FVMProblem from a given CellProblem.\n\n\n\n\n\n","category":"type"},{"location":"#MovingBoundaryProblems1D.MBProblem","page":"Home","title":"MovingBoundaryProblems1D.MBProblem","text":"MBProblem(prob::CellProblem,\n    mesh_points=copy(prob.initial_condition);\n    diffusion_function=:continuum,\n    diffusion_parameters=nothing,\n    reaction_function=:continuum,\n    reaction_parameters=nothing,\n    moving_boundary_function=:continuum,\n    moving_boundary_parameters=nothing,\n    rhs_function = :continuum,\n    rhs_parameters=nothing,\n    proliferation=false)\n\nConstructs an MBProblem from a given CellProblem. \n\n\n\n\n\n","category":"type"},{"location":"#EpithelialDynamics1D.cell_densities-Union{Tuple{AbstractVector{T}}, Tuple{T}} where T<:Number","page":"Home","title":"EpithelialDynamics1D.cell_densities","text":"cell_densities(cell_positions::AbstractVector{T}) where {T<:Number}\n\nCompute the cell densities from the cell positions. The returned vector q has  length(q) = length(cell_positions) - 1, with q[i] the density of the ith cell (cell_positions[i], cell_positions[i+1]) given by  1/(cell_positions[i+1] - cell_positions[i]).\n\nIf you want estimates of the densities at each node rather than for a cell,  see node_densities.\n\n\n\n\n\n","category":"method"},{"location":"#EpithelialDynamics1D.cell_midpoints-Union{Tuple{AbstractVector{T}}, Tuple{T}} where T<:Number","page":"Home","title":"EpithelialDynamics1D.cell_midpoints","text":"cell_midpoints(cell_positions::AbstractVector{T}) where {T<:Number}\n\nCompute the midpoints of the cells from the cell positions. The returned vector x has length(x) = length(cell_positions) - 1, with x[i] the midpoint of the ith cell (cell_positions[i], cell_positions[i+1]) given by  0.5(cell_positions[i] + cell_positions[i+1]).\n\n\n\n\n\n","category":"method"},{"location":"#EpithelialDynamics1D.cell_numbers-Tuple{SciMLBase.EnsembleSolution}","page":"Home","title":"EpithelialDynamics1D.cell_numbers","text":"cell_numbers(sol::EnsembleSolution; indices = eachindex(sol), alpha=0.05)\n\nComputes summary statistics for the cell numbers from an EnsembleSolution to a CellProblem.\n\nArguments\n\nsol::EnsembleSolution: The ensemble solution to a CellProblem.\n\nKeyword Arguments\n\nindices = eachindex(sol): The indices of the cell simulations to consider.\nalpha::Float64 = 0.05: The significance level for the confidence intervals.\n\nOutputs\n\nN::Vector{Vector{Int}}: The cell numbers for each cell simulation.\nmeans::Vector{Float64}: The mean cell numbers for each cell simulation.\nlowers::Vector{Float64}: The lower bounds of the confidence intervals for the cell numbers for each cell simulation.\nuppers::Vector{Float64}: The upper bounds of the confidence intervals for the cell numbers for each cell simulation.\n\n\n\n\n\n","category":"method"},{"location":"#EpithelialDynamics1D.continuum_limit","page":"Home","title":"EpithelialDynamics1D.continuum_limit","text":"continuum_limit(prob::CellProblem,\n    mesh_points=copy(prob.initial_condition);\n    proliferation=false)\n\nReturns the FVMProblem (if the right boundary is fixed) or  the MBProblem (if the right boundary is free) that defines the  continuum limit to the CellProblem prob. The argument  mesh_points defines the mesh points to use for the PDE problem,  which can be a grid of points or a number of points to use. If the problem  should be considered to have proliferation, set proliferation = true.\n\n\n\n\n\n","category":"function"},{"location":"#EpithelialDynamics1D.get_knots","page":"Home","title":"EpithelialDynamics1D.get_knots","text":"get_knots(sol, num_knots = 500; indices = eachindex(sol))\n\nComputes knots for each time, covering the extremum of the cell positions across all  cell simulations. You can restrict the simultaions to consider using the indices.\n\n\n\n\n\n","category":"function"},{"location":"#EpithelialDynamics1D.integrate_pde-Tuple{Any, Any}","page":"Home","title":"EpithelialDynamics1D.integrate_pde","text":"integrate_pde(sol, i)\n\nIntegrate the PDE solution sol at time i, giving an approximation of the number of cells at sol.t[i].\n\n\n\n\n\n","category":"method"},{"location":"#EpithelialDynamics1D.leading_edges-Tuple{SciMLBase.EnsembleSolution}","page":"Home","title":"EpithelialDynamics1D.leading_edges","text":"leading_edges(sol::EnsembleSolution; indices = eachindex(sol), alpha=0.05)\n\nComputes summary statistics for the leading edges from an EnsembleSolution to a CellProblem.\n\nArguments\n\nsol::EnsembleSolution: The ensemble solution to a CellProblem.\n\nKeyword Arguments\n\nindices = eachindex(sol): The indices of the cell simulations to consider.\nalpha::Float64 = 0.05: The significance level for the confidence intervals.\n\nOutputs\n\nL::Vector{Vector{Float64}}: The leading edges for each cell simulation.\nmeans::Vector{Float64}: The mean leading edges for each cell simulation.\nlowers::Vector{Float64}: The lower bounds of the confidence intervals for the leading edges for each cell simulation.\nuppers::Vector{Float64}: The upper bounds of the confidence intervals for the leading edges for each cell simulation.\n\n\n\n\n\n","category":"method"},{"location":"#EpithelialDynamics1D.node_densities-Tuple{SciMLBase.EnsembleSolution}","page":"Home","title":"EpithelialDynamics1D.node_densities","text":"node_densities(sol::EnsembleSolution; num_knots=500, knots=get_knots(sol, num_knots), alpha=0.05, interp_fnc=(u, t) -> LinearInterpolation{true}(u, t))\n\nComputes summary statistics for the node densities from an EnsembleSolution to a CellProblem.\n\nArguments\n\nsol::EnsembleSolution: The ensemble solution to a CellProblem.\n\nKeyword Arguments\n\nindices = eachindex(sol): The indices of the cell simulations to consider. \nnum_knots::Int = 500: The number of knots to use for the spline interpolation.\nknots::Vector{Vector{Float64}} = get_knots(sol, num_knots; indices): The knots to use for the spline interpolation.\nalpha::Float64 = 0.05: The significance level for the confidence intervals.\ninterp_fnc = (u, t) -> LinearInterpolation{true}(u, t): The function to use for constructing the interpolant.\n\nOutputs\n\nq::Vector{Vector{Vector{Float64}}}: The node densities for each cell simulation.\nr::Vector{Vector{Vector{Float64}}}: The cell positions for each cell simulation.\nmeans::Vector{Vector{Float64}}: The mean node densities for each cell simulation.\nlowers::Vector{Vector{Float64}}: The lower bounds of the confidence intervals for the node densities for each cell simulation.\nuppers::Vector{Vector{Float64}}: The upper bounds of the confidence intervals for the node densities for each cell simulation.\nknots::Vector{Vector{Float64}}: The knots used for the spline interpolation.\n\n\n\n\n\n","category":"method"},{"location":"#EpithelialDynamics1D.node_densities-Union{Tuple{AbstractVector{T}}, Tuple{T}} where T<:Number","page":"Home","title":"EpithelialDynamics1D.node_densities","text":"node_densities(cell_positions::AbstractVector{T}) where {T<:Number}\n\nCompute the cell densities from the cell positions, assigning a density to each node.  The ith density is given by 2/(cell_positions[i+1] - cell_positions[i-1]) if i is not the first or last node, 1/(cell_positions[i+1] - cell_positions[i]) if i is the first node, and 1/(cell_positions[i] - cell_positions[i-1]) if i is the last node.\n\nIf you want estimates of the densities for each cell rather than at each node, see cell_densities.\n\n\n\n\n\n","category":"method"}]
}
