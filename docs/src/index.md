```@meta
CurrentModule = EpithelialDynamics1D
```

# EpithelialDynamics1D

Documentation for [EpithelialDynamics1D](https://github.com/DanielVandH/EpithelialDynamics1D.jl).

This is a package for simulating epithelial dynamics in one dimension, supporting cell migration and cell proliferation, implementing the model in [Baker et al. (2019)](https://doi.org/10.1016/j.jtbi.2018.12.025). In this model, cells are represented as intervals between points, with the $i$th cell given by $(x_i, x_{i+1})$, $i=1,\ldots,n-1$. The node positions $x_i$ are governed by the differential equation

```math
\eta\dfrac{\mathrm dx_i}{\mathrm dt} = F\left(x_i - x_{i-1}\right) - F\left(x_{i+1} - x_i\right), \quad i=2,\ldots,n-1,
```

assuming $x_1 < x_2 < \cdots < x_n$. If the left boundary is fixed, then $\mathrm dx_1/\mathrm dt = 0$, otherwise $\eta\mathrm dx_1/\mathrm dt = -F(x_2 - x_1)$. Similarly, if the right boundary is fixed then $\mathrm dx_n/\mathrm dt = 0$, otherwise $\eta\mathrm dx_n/\mathrm dt = F(x_n - x_{n-1})$. The parameter $\eta$ is called the _damping constant_ or the _viscosity coefficient_ and $F$ is the _force law_. In this model, cells are modelled as springs, whose forces are governed by this force law $F$. 

The proliferation mechanism that we use assumes that at most one cell can divide at a time. Let $C_i(t)$ denote the event that the $i$th cell divides in the time interval $[t, t+\mathrm dt)$. We assume that $\Pr[C_i(t)] = G_i\,\mathrm dt$, where $G_i = G(|x_{i+1} - x_i|)$ for some proliferation law $G$. When this division occurs we place a new position at $(x_i + x_{i+1})/2$, adjusting the indices of the $x_i$ accordingly so that they remain sorted. We implement this mechanism by attempting a proliferation event every $\mathrm dt = \Delta t$ units of time, so that the probability of a proliferation event occuring at the time $t$ is $\sum_{i=1}^{n-1} G_i\Delta t$, and the probability that the $i$th cell proliferates, given that a proliferation event does occur, is $G_i/\sum_{j=1}^{n-1} G_j$. [DiffEqCallbacks.jl](https://github.com/SciML/DiffEqCallbacks.jl) is used to implement the callback that attempts this event periodically while simulating the problem.

Examples of how to simulate these systems are given in the sidebar.

```@index
```

```@autodocs
Modules = [EpithelialDynamics1D]
```
