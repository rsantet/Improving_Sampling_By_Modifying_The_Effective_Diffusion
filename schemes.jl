include("input.jl")

function hamiltonian(q, p, ∇z, unit_vector, distance, α, κ_α)

    z = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)
    a_z = a(z, α)

    return V(q, unit_vector, distance) - log(a_z) / 2 + κ_α * dot(p, p) / 2 + κ_α * (a_z - 1) / 2 * dot(reduce(vcat, ∇z), reduce(vcat, p[1:2]))^2 / square_norm_∇z
end

function ∇p_hamiltonian(q, p, ∇z, unit_vector, distance, α, κ_α)
    z = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)
    res = κ_α * p
    res[1:2] .+= κ_α * (a(z, α) - 1) * dot(∇z, p[1:2]) * ∇z / square_norm_∇z
    return res
end

function ∇q_hamiltonian(q, p, F, ∇z, ∇2z, unit_vector, distance, α, κ_α)
    compute_forces!(q, F, unit_vector, distance)
    res = -F
    z = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)
    ∇2ξ!(q, ∇2z, unit_vector, distance)
    a_z = a(z, α)
    ∇a_z = ∇a(z, α)
    dot_z = dot(reduce(vcat, ∇z), reduce(vcat, p[1:2]))
    vec_z = ∇2z * p[1:2]
    res[1:2] .+= (
        -1 / 2 * ∇a_z / a_z * ∇z + κ_α / 2 * dot_z^2 * ∇a_z * ∇z / square_norm_∇z + κ_α * (a_z - 1) * dot_z * vec_z / square_norm_∇z
    )
    return res
end

function ∇2qp_hamiltonian(q, p, ∇z, ∇2z, unit_vector, distance, α, κ_α)

    z = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)
    ∇2ξ!(q, ∇2z, unit_vector, distance)
    a_z = a(z, α)
    ∇a_z = ∇a(z, α)
    dot_z = dot(∇z, p[1:2])
    vec_z = ∇2z * p[1:2]

    return κ_α * (
        dot_z * ∇a_z * outer_product(∇z, ∇z) + (a_z - 1) * outer_product(∇z, vec_z) + (a_z - 1) * dot_z * ∇2z
    ) / square_norm_∇z
end

function GSV_momenta_update!(q, p, F, unit_vector, distance, ∇z, ∇2z, α, κ_α, Δt)

    # Forces should already be computed at this point
    z = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)
    ∇2ξ!(q, ∇2z, unit_vector, distance)
    a_z = a(z, α)
    ∇a_z = ∇a(z, α)
    outer_∇z = outer_product(∇z, ∇z)
    norm_init = norm(p[1:2])

    # Initialize the Newton sequence
    old_p = deepcopy(p)
    dot_z = dot(∇z, old_p[1:2])
    vec_z = ∇2z * old_p[1:2]
    new_p = old_p + Δt / 2 * F # Exact execpt for first four components
    new_p[1:2] -= Δt / 2 * (
        (
            -∇a_z / a_z / 2 + κ_α * dot_z^2 * ∇a_z / square_norm_∇z / 2
        ) * ∇z
        +
        κ_α * (a_z - 1) * dot_z * vec_z / square_norm_∇z
    )

    # Compute value of the function
    dot_z = dot(∇z, new_p[1:2])
    vec_z = ∇2z * new_p[1:2]
    g_new_p = new_p - p - Δt / 2 * F
    g_new_p[1:2] .+= Δt / 2 * (
        (
            -∇a_z / a_z / 2 + κ_α / 2 * ∇a_z * dot_z^2 / square_norm_∇z
        ) * ∇z
        +
        κ_α * (a_z - 1) * dot_z * vec_z / square_norm_∇z
    )

    # println("INITIAL VALUE")
    # display(g_new_p)
    # println()

    # Initialize counter
    counter = 0

    while counter < max_iterations_newton && (
        norm(g_new_p[1:2]) > tol_root_newton ||
        norm(new_p[1:2] - old_p[1:2]) > tol_cauchy_newton * norm_init
    )

        # Build Jacobian of g
        ∇g = I + κ_α * Δt / 2 * (
            dot_z * ∇a_z * outer_∇z + (a_z - 1) * outer_product(vec_z, ∇z) + (a_z - 1) * dot_z * ∇2z
        ) / square_norm_∇z

        rank_∇g = try
            rank(∇g)
        catch
            return false
        end

        if rank_∇g < 4
            return false
        end

        # Solve linear system
        res = ∇g \ g_new_p[1:2]

        # Update the sequence
        old_p = deepcopy(new_p)
        # new_p -= g_new_p
        new_p[1:2] = old_p[1:2] - res

        # Update the value of g
        dot_z = dot(∇z, new_p[1:2])
        vec_z = ∇2z * new_p[1:2]
        g_new_p = new_p - p - Δt / 2 * F
        g_new_p[1:2] .+= Δt / 2 * (
            (
                -∇a_z / a_z / 2 + κ_α * ∇a_z / 2 * dot_z^2 / square_norm_∇z
            ) * ∇z
            +
            κ_α * (a_z - 1) * dot_z * vec_z / square_norm_∇z
        )

        # Increase counter
        counter += 1
    end

    # If the method has converged
    if norm(g_new_p[1:2]) < tol_root_newton && norm(new_p[1:2] - old_p[1:2]) < tol_cauchy_newton * norm_init
        p .= new_p

        return true
    end

    # println("Not converged")
    # println("ENDING VALUE")
    # display(g_new_p)
    # println()
    # println(norm(g_new_p), " ", tol_root_newton)
    # println(norm(new_p-old_p), " ", tol_cauchy_newton*norm_init)
    # println()

    return false
end

function GSV_position_update!(q, p, unit_vector, distance, ∇z, ∇2z, α, κ_α, Δt)

    z_init = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)
    ∇2ξ!(q, ∇2z, unit_vector, distance)
    ∇z_init = deepcopy(∇z)
    dot_z_init = dot(∇z_init, p[1:2])
    a_z_init = a(z_init, α)
    norm_init = norm(q[1:2])

    # Initialize the Newton sequence
    old_q = deepcopy(q)
    new_q = old_q + Δt * κ_α * p # Exact for all components except the first four
    new_q[1:2] .+= κ_α * Δt * (a_z_init - 1) * dot_z_init * ∇z_init / square_norm_∇z

    # Compute value of the function
    z = ξ(new_q, unit_vector, distance)
    a_z = a(z, α)
    ∇a_z = ∇a(z, α)
    ∇ξ!(new_q, ∇z, unit_vector, distance)
    ∇2ξ!(new_q, ∇2z, unit_vector, distance)
    dot_z = dot(∇z, p[1:2])
    h_new_q = new_q - q - κ_α * Δt * p
    h_new_q[1:2] -= Δt * κ_α / 2 / square_norm_∇z * (
        (a_z_init - 1) * dot_z_init * ∇z_init + (a_z - 1) * dot_z * ∇z
    )

    # Initialize counter
    counter = 0

    while counter < max_iterations_newton && (norm(h_new_q[1:2]) > tol_root_newton || norm(new_q[1:2] - old_q[1:2]) > tol_cauchy_newton * norm_init)

        # Build Jacobian matrix of h
        ∇h = I - κ_α * Δt / 2 / square_norm_∇z * (
            dot_z * ∇a_z * outer_product(∇z, ∇z) + (a_z - 1) * outer_product(∇z, ∇2z * p[1:2]) + (a_z - 1) * dot_z * ∇2z
        )

        rank_∇h = try
            rank(∇h)
        catch
            return false
        end

        if rank_∇h < 4
            return false
        end

        # Solve linear system ∇h(q^i)(q^{i+1}-q^i)=-h(q^i)
        res = ∇h \ h_new_q[1:2] # (2,)

        # Update the sequence
        old_q = deepcopy(new_q)
        # new_q -= h_new_q
        new_q[1:2] = old_q[1:2] - res

        # Update the value of h
        z = ξ(new_q, unit_vector, distance)
        ∇ξ!(new_q, ∇z, unit_vector, distance)
        ∇2ξ!(new_q, ∇2z, unit_vector, distance)
        a_z = a(z, α)
        ∇a_z = ∇a(z, α)
        dot_z = dot(∇z, p[1:2])
        h_new_q = new_q - q - κ_α * Δt * p
        h_new_q[1:2] -= Δt * κ_α / 2 / square_norm_∇z * (
            (a_z_init - 1) * dot_z_init * ∇z_init + (a_z - 1) * dot_z * ∇z
        )

        counter += 1

    end

    # If the method has converged
    if norm(h_new_q[1:2]) < tol_root_newton && norm(new_q[1:2] - old_q[1:2]) < tol_cauchy_newton * norm_init

        # println("Newton convergence after $counter")
        q .= new_q
        return true
    end

    return false

end

function GSV!(q, p, F, unit_vector, distance, ∇z, ∇2z, α, κ_α, Δt, rejection_counter::Counter, mode::String)

    q_init = deepcopy(q)
    p_init = deepcopy(p)
    F_init = deepcopy(F)

    # Implicit over momenta
    convergence = GSV_momenta_update!(q, p, F, unit_vector, distance, ∇z, ∇2z, α, κ_α, Δt)


    if !convergence
        if mode == "Forward"
            rejection_counter.forward_no_convergence_momenta_update += 1
        elseif mode == "Backward"
            rejection_counter.backward_no_convergence_momenta_update += 1
        end
        q .= q_init
        p .= p_init
        F .= F_init
        return false
    end

    # Implicit over position
    convergence = GSV_position_update!(q, p, unit_vector, distance, ∇z, ∇2z, α, κ_α, Δt)


    if !convergence
        if mode == "Forward"
            rejection_counter.forward_no_convergence_position_update += 1
        elseif mode == "Backward"
            rejection_counter.backward_no_convergence_position_update += 1
        end
        q .= q_init
        p .= p_init
        F .= F_init
        return false
    end

    # Explicit over momenta
    compute_forces!(q, F, unit_vector, distance) # F now holds the forces of the proposal
    z = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)
    ∇2ξ!(q, ∇2z, unit_vector, distance)
    a_z = a(z, α)
    ∇a_z = ∇a(z, α)
    dot_z = dot(∇z, p[1:2])
    vec_z = ∇2z * p[1:2]

    p .+= Δt / 2 * F
    p[1:2] .-= Δt / 2 * (
        (
            -∇a_z / a_z / 2 + κ_α * dot_z^2 * ∇a_z / square_norm_∇z / 2
        ) * ∇z
        +
        κ_α * (a_z - 1) * dot_z * vec_z / square_norm_∇z
    )


    # z = ξ(q, unit_vector, distance) # Check if CV value exceeds the range [z_min,z_max]

    return true

end

function one_step_RMHMC!(q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, inv_sqrt_D, α, κ_α, Δt, rejection_counter::Counter)

    # Resample momenta
    z = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)
    projection_matrix .= outer_product(∇z, ∇z) / square_norm_∇z
    inv_sqrt_D .= (I + (1 / sqrt(a(z, α)) - 1) * projection_matrix) / sqrt(κ_α)
    p = randn(SVector{2,Float64}, n_particles) / sqrt(β)
    p[1:2] .= inv_sqrt_D * p[1:2]


    q_init = deepcopy(q)
    F_init = deepcopy(F)
    H_init = hamiltonian(q, p, ∇z, unit_vector, distance, α, κ_α)

    # Integrate the Hamiltonian dynamics
    convergence = GSV!(q, p, F, unit_vector, distance, ∇z, ∇2z, α, κ_α, Δt, rejection_counter, "Forward")

    if !convergence
        q .= q_init
        F .= F_init
        return
    end

    rev_q = deepcopy(q)
    rev_p = -deepcopy(p)
    rev_F = deepcopy(F) # Holds the forces of the proposal, no need to recompute


    convergence = GSV!(rev_q, rev_p, rev_F, unit_vector, distance, ∇z, ∇2z, α, κ_α, Δt, rejection_counter, "Backward")

    if !convergence
        q .= q_init
        F .= F_init
        return
    end

    rev_norm = norm((rev_q-q_init)[1:2])
    rev_check = tol_reversibility * norm(q_init[1:2])
    rev = rev_norm < rev_check

    if !rev
        # println("NO REV")
        # println(rev_norm)
        # println(rev_check)
        # println()
        rejection_counter.no_reversibility += 1
        q .= q_init
        F .= F_init
        return
    end

    H = hamiltonian(q, p, ∇z, unit_vector, distance, α, κ_α)

    if log(rand()) > β * (H_init - H)
        rejection_counter.MH_rejection += 1
        q .= q_init
        F .= F_init
        return
    end

end

function OU_update!(q, p, ∇z, unit_vector, distance, projection_matrix, γ, Δt, α, κ_α)

    # Compute the projection matrix
    z = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)
    projection_matrix .= outer_product(∇z, ∇z) / square_norm_∇z
    a_z = a(z, α)
    tmp = 1 / (1 + Δt / 4 * γ * κ_α)

    # Update of the first two components
    p[1:2] .= (
        I * tmp + (1 / (1 + Δt / 4 * γ * κ_α * a_z) - tmp) * projection_matrix
    ) * (
        (1 - Δt / 4 * γ * κ_α) * p[1:2] + sqrt(γ * Δt / β) * randn(SVector{2,Float64}, 2) - Δt / 4 * γ * κ_α * (a_z - 1) * projection_matrix * p[1:2]
    )

    # Update of the remaining components
    p[3:end] .= (
        (1 - Δt / 4 * γ * κ_α) * p[3:end] + sqrt(γ * Δt / β) * randn(SVector{2,Float64}, n_particles - 2)
    ) * tmp

end

function one_step_RMGHMC!(q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, α, κ_α, Δt, γ, rejection_counter::Counter)

    # OU update with midpoint Euler scheme
    OU_update!(q, p, ∇z, unit_vector, distance, projection_matrix, γ, Δt, α, κ_α)

    q_init = deepcopy(q)
    p_init = deepcopy(p)
    F_init = deepcopy(F)
    H_init = hamiltonian(q, p, ∇z, unit_vector, distance, α, κ_α)

    # Integrate the Hamiltonian dynamics
    forward_convergence = GSV!(q, p, F, unit_vector, distance, ∇z, ∇2z, α, κ_α, Δt, rejection_counter, "Forward")

    if forward_convergence

        rev_q = deepcopy(q)
        rev_p = -deepcopy(p)
        rev_F = deepcopy(F) # Holds the forces of the proposal, no need to recompute

        backward_convergence = GSV!(rev_q, rev_p, rev_F, unit_vector, distance, ∇z, ∇2z, α, κ_α, Δt, rejection_counter, "Backward")

        if backward_convergence

            rev_norm = norm((rev_q-q_init)[1:2])
            rev_check = tol_reversibility * norm(q_init[1:2])
            rev = rev_norm < rev_check

            if rev

                H = hamiltonian(q, p, ∇z, unit_vector, distance, α, κ_α)

                if log(rand()) > β * (H_init - H)
                    rejection_counter.MH_rejection += 1
                    q .= q_init
                    p .= -p_init
                    F .= F_init
                end

            else # No reversibility
                rejection_counter.no_reversibility += 1
                q .= q_init
                p .= -p_init
                F .= F_init

            end
        else # No backward convergence
            q .= q_init
            p .= -p_init
            F .= F_init
        end
    else # No forward convergence

        q .= q_init
        p .= -p_init
        F .= F_init

    end

    # OU update with midpoint Euler scheme
    OU_update!(q, p, ∇z, unit_vector, distance, projection_matrix, γ, Δt, α, κ_α)

end

function one_step_MALA!(q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, div_D, inv_D, Δt, α, κ_α, rejection_counter::Counter)

    # Copies in case the move is rejected
    old_q = deepcopy(q)
    old_F = deepcopy(F)
    old_sqrt_D = deepcopy(sqrt_D)
    old_D = deepcopy(D)
    old_div_D = deepcopy(div_D)

    old_z = ξ(old_q, unit_vector, distance)
    old_a_z = a(old_z, α)
    old_V = V(q, unit_vector, distance)


    noise = randn(SVector{2,Float64}, n_particles)

    # First two components
    q[1:2] .+= (
        (old_D * F[1:2] + old_div_D / β) * Δt
        +
        sqrt(2 * Δt / β) * old_sqrt_D * noise[1:2]
    )
    # Remaining components
    q[3:end] .+= (Δt * κ_α * F[3:end] + sqrt(2 * Δt * κ_α / β) * noise[3:end])


    # Compute MH ratio
    z = ξ(q, unit_vector, distance)
    a_z = a(z, α)

    compute_sqrt_div_inv_diffusion!(q, ∇z, projection_matrix, D, sqrt_D, div_D, inv_D, unit_vector, distance, α, κ_α) # Now holds the values for the proposal

    new_V = V(q, unit_vector, distance)
    compute_forces!(q, F, unit_vector, distance) # F now stores the forces of the proposal

    # First four components
    drift_first_four_components = (D * F[1:2] + div_D / β) * Δt
    G_proposal_1 = dot(old_q[1:2] - q[1:2] - drift_first_four_components, inv_D * (old_q[1:2] - q[1:2] - drift_first_four_components))

    # Remaining components
    G_proposal_2 = dot(old_q[3:end] - q[3:end] - F[3:end] * Δt * κ_α, old_q[3:end] - q[3:end] - F[3:end] * Δt * κ_α) / κ_α

    G_proposal = (G_proposal_1 + G_proposal_2) * β / 4 / Δt

    log_r = β * (old_V - new_V) + log(old_a_z / a_z) / 2 + square_norm(noise) / 2 - G_proposal

    if log(rand()) > log_r # Rejection
        q .= old_q
        F .= old_F
        sqrt_D .= old_sqrt_D
        D .= old_D
        div_D .= old_div_D

        rejection_counter.MH_rejection += 1
        reaction_coordinate_value = old_z
    else
        reaction_coordinate_value = z
    end

    return reaction_coordinate_value
end

function one_step_adaptive_MALA!(q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, div_D, inv_D, Δt, α, κ_α, sqrt_κ_α, MeanForceUpdate::Matrix, MeanForce::Matrix, FreeEnergy::Vector, rejection_counter::Counter)

    # Copies in case the move is rejected
    old_q = deepcopy(q)
    old_F = deepcopy(F)
    old_sqrt_D = deepcopy(sqrt_D)
    old_D = deepcopy(D)
    old_div_D = deepcopy(div_D)

    old_V = V(q, unit_vector, distance)
    old_z = ξ(old_q, unit_vector, distance)
    old_index = get_index_RC(old_z)

    old_MF = get_mean_force(old_index, MeanForce)
    old_FE = get_free_energy(old_index, old_z, FreeEnergy, MeanForce)
    old_a_z = adaptive_a(old_FE, α)

    noise = randn(SVector{2,Float64}, n_particles)

    # First two components
    q[1:2] .+= (
        (old_D * F[1:2] + old_div_D / β) * Δt
        +
        sqrt(2 * Δt / β) * old_sqrt_D * noise[1:2]
    )
    # Remaining components
    q[3:end] .+= (Δt * κ_α * F[3:end] + sqrt(2 * Δt * κ_α / β) * noise[3:end])


    # Compute MH ratio
    z = ξ(q, unit_vector, distance)
    index = get_index_RC(z)
    MF = get_mean_force(index, MeanForce)
    FE = get_free_energy(index, z, FreeEnergy, MeanForce)
    a_z = adaptive_a(FE, α)

    adaptive_compute_sqrt_div_inv_diffusion!(q, ∇z, projection_matrix, D, sqrt_D, div_D, inv_D, unit_vector, distance, MF, FE, α, κ_α, sqrt_κ_α) # Now holds the values for the proposal

    new_V = V(q, unit_vector, distance)
    compute_forces!(q, F, unit_vector, distance) # F now stores the forces of the proposal

    # First four components
    drift_first_four_components = (D * F[1:2] + div_D / β) * Δt
    G_proposal_1 = dot(old_q[1:2] - q[1:2] - drift_first_four_components, inv_D * (old_q[1:2] - q[1:2] - drift_first_four_components))

    # Remaining components
    G_proposal_2 = dot(old_q[3:end] - q[3:end] - F[3:end] * Δt * κ_α, old_q[3:end] - q[3:end] - F[3:end] * Δt * κ_α) / κ_α

    G_proposal = (G_proposal_1 + G_proposal_2) * β / 4 / Δt

    log_r = β * (old_V - new_V) + log(old_a_z / a_z) / 2 + square_norm(noise) / 2 - G_proposal

    if log(rand()) > log_r # Rejection
        q .= old_q
        F .= old_F
        sqrt_D .= old_sqrt_D
        D .= old_D
        div_D .= old_div_D

        rejection_counter.MH_rejection += 1
        reaction_coordinate_value = old_z
        MF = old_MF
        FE = old_FE
        index = old_index
        new_V = old_V
    else
        reaction_coordinate_value = z
    end

    # Update the mean force update
    # Compute the local mean force
    length_PBC!(q, unit_vector, distance, 1, 2)
    lmf = (dot(F[2] - F[1], unit_vector[1]) - 2 / β / distance[1]) * w
    if index != 0
        MeanForceUpdate[index, 1] += 1
        MeanForceUpdate[index, 2] += lmf
    end

    return reaction_coordinate_value
end