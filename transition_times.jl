include("input.jl")
include("schemes.jl")

function compute_transition_times(scheme, Δt, n_transitions::Int, set_tolerance, α; n_loggers=nothing, γ=1.0)

    data_dir = "data/$scheme/transition_times/$(round(α, sigdigits=5))/$(round(Δt, sigdigits=5))/"

    mkpath(data_dir)
    tt_fp = data_dir * "tt.jld2"
    n_iterations_fp = data_dir * "n_iterations.jld2"

    if isfile(tt_fp) && isfile(n_iterations_fp)
        n_iterations_vec = load_object(n_iterations_fp)
        tt_vec = load_object(tt_fp)

        len_n_iterations = length(n_iterations_vec)
        len_tt = length(tt_vec)

        if len_n_iterations >= n_transitions && len_tt >= n_transitions
            println("Found $len_n_iterations/$len_tt transitions")
            return tt_vec[1:n_transitions], n_iterations_vec[1:n_transitions]
        end
    end

    println("Running for Δt = ", Δt, " and α = ", α)

    ## Run the long trajectory
    # Initial condition
    initial_z = 0.0

    if scheme == "MALA"
        q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, inv_D, div_D, κ_α = initialize_OL(initial_z, α)
    elseif scheme == "adaptive_MALA"
        q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, inv_D, div_D, MeanForce, FreeEnergy, MeanForceUpdate, κ_α, sqrt_κ_α = initialize_adaptive_OL(initial_z, α)
    elseif scheme == "RMHMC" || scheme == "RMGHMC"
        q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, inv_sqrt_D, κ_α = initialize_HMC(initial_z, α)
    end

    rejection_counter = Counter(0, 0, 0, 0, 0, 0)

    compact_form = true # Initial configuration is in compact form
    time = 0.0 # Keeping track of the time passed
    n_iterations = 0 # Overall iterations
    n_iterations_counter = 0 # Keeping track of the number of iterations in-between transitions
    n_transitions_counter = 0 # Keeping track of the number of transitions

    tt_vec = zeros(n_transitions)
    n_iterations_vec = zeros(Int, n_transitions)

    while n_transitions_counter < n_transitions

        time += Δt
        n_iterations_counter += 1
        n_iterations += 1

        if !isnothing(n_loggers)
            if n_iterations_counter % n_loggers == 0
                z = ξ(q, unit_vector, distance)
                println("$scheme, α = ", round(α, sigdigits=5), "\t Δt = ", round(Δt, sigdigits=5), "\t Transition\t", n_transitions_counter, "\tIteration\t", n_iterations_counter, "\t RC:\t", z)
            end
        end

        # Update
        if scheme == "RMHMC"
            one_step_RMHMC!(q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, inv_sqrt_D, α, κ_α, Δt, rejection_counter)
        elseif scheme == "RMGHMC"
            one_step_RMGHMC!(q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, α, κ_α, Δt, γ, rejection_counter)
        elseif scheme == "MALA"
            one_step_MALA!(q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, div_D, inv_D, Δt, α, κ_α, rejection_counter)
        elseif scheme == "adaptive_MALA"

            if n_iterations % n_mean_force_update == 0
                update_mean_force!(MeanForceUpdate, MeanForce)
                update_free_energy!(MeanForce, FreeEnergy)
                κ_α, sqrt_κ_α = update_κ_α(FreeEnergy, α)
            end

            one_step_adaptive_MALA!(q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, div_D, inv_D, Δt, α, κ_α, sqrt_κ_α, MeanForce, MeanForce, FreeEnergy, rejection_counter)
        end

        # Check if a transition has occurred
        if compact_form
            if is_in_expanded_form(q, unit_vector, distance, set_tolerance)

                compact_form = false
                n_transitions_counter += 1


                println("$scheme, α = ", round(α, sigdigits=5), "\t Δt = ", round(Δt, sigdigits=5), "\t Transition 0 -> 1 has occurred after ", n_iterations_counter, " iterations. Remaining transitions: ", n_transitions - n_transitions_counter)

                tt_vec[n_transitions_counter] = time
                n_iterations_vec[n_transitions_counter] = n_iterations_counter
                time = 0.0
                n_iterations_counter = 0
            end
        else
            if is_in_compact_form(q, unit_vector, distance, set_tolerance)

                compact_form = true
                n_transitions_counter += 1

                println("$scheme, α = ", round(α, sigdigits=5), "\t Δt = ", round(Δt, sigdigits=5), "\t Transition 1 -> 0 has occurred after ", n_iterations_counter, " iterations. Remaining transitions: ", n_transitions - n_transitions_counter)

                tt_vec[n_transitions_counter] = time
                n_iterations_vec[n_transitions_counter] = n_iterations_counter
                time = 0.0
                n_iterations_counter = 0
            end
        end
    end

    save_object(tt_fp, tt_vec)
    save_object(n_iterations_fp, n_iterations_vec)
    return tt_vec, n_iterations_vec
end

# n_transitions = 100
# Δt = 1e-3
# n_loggers = floor(Int, 1 / Δt)
# α = 0.0
# γ = 1.0
# scheme = "adaptive_MALA"

# res = compute_transition_times(scheme, Δt, n_transitions, set_tolerance, α; γ=γ)
# using Statistics
# println(mean(res[2]))