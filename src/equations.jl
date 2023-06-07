function cell_odes!(dr, r, prob, t)
    η = prob.damping_constant
    η⁻¹ = inv(η)
    fix_left = prob.fix_left
    fix_right = prob.fix_right
    F = prob.force_law
    p = prob.force_law_parameters
    for i in (firstindex(dr)+1):(lastindex(dr)-1)
        δᵢ₋₁ = r[i] - r[i-1]
        δᵢ = r[i+1] - r[i]
        Fᵢ₋₁ = F(δᵢ₋₁, p)
        Fᵢ = F(δᵢ, p)
        dr[i] = η⁻¹ * (Fᵢ₋₁ - Fᵢ)
    end
    if !fix_left
        δ₁ = r[begin+1] - r[begin]
        F₁ = F(δ₁, p)
        dr[begin] = -η⁻¹ * F₁
    else
        dr[begin] = zero(dr[begin])
    end
    if !fix_right
        δₙ₋₁ = r[end] - r[end-1]
        Fₙ₋₁ = F(δₙ₋₁, p)
        dr[end] = η⁻¹ * Fₙ₋₁
    else
        dr[end] = zero(dr[end])
    end
    return nothing
end

function get_nnz_itr(i, n)
    if i == 1
        return (1, 2)
    elseif 1 < i < n
        return (i - 1, i, i + 1)
    else # if i == n 
        return (n - 1, n)
    end
end
jacobian_sparsity(prob::CellProblem) = jacobian_sparsity(prob.initial_condition)
function jacobian_sparsity(pts)
    n = length(pts)
    num_nnz = 3n - 2
    I = zeros(Int, num_nnz)
    J = zeros(Int, num_nnz)
    V = ones(num_nnz)
    ctr = 1
    for i in 1:n
        itr = get_nnz_itr(i, n)
        for j in itr
            I[ctr] = i
            J[ctr] = j
            ctr += 1
        end
    end
    return sparse(I, J, V)
end