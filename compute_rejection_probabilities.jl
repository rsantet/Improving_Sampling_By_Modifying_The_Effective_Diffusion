using JLD2
include("input.jl")
include("schemes.jl")

using Base.Threads

function rejection_probabilities_Δt(scheme, n_iterations, Δt, α; γ=1.0)

    α_str = round(α, sigdigits=5)
    Δt_str = round(Δt, sigdigits=5)
    dir_path = "data/$scheme/rejection_probabilities/$n_iterations/$α_str/$Δt_str/"
    if isfile(dir_path * "proba_global_rejection.jld2")
        return
    end

    # Initialize the system
    initial_z = 0.0
    q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, inv_sqrt_D, κ_α = initialize_HMC(initial_z, α)

    # thermalizing run

    println("$α_str\t $Δt_str\t thermalizing run")
    rejection_counter_tr = Counter(0, 0, 0, 0, 0, 0)



    for i in floor(Int, 1 / Δt)

        if i % 100000 == 0
            println("$α_str\t $Δt_str\t  $i / $n_iterations")
            flush(stdout)
        end

        if scheme == "RMHMC"
            one_step_RMHMC!(q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, inv_sqrt_D, α, κ_α, Δt, rejection_counter_tr)
        elseif scheme == "RMGHMC"
            one_step_RMGHMC!(q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, α, κ_α, Δt, γ, rejection_counter_tr)
        end

    end


    println("$α_str\t $Δt_str\t Computing rejection probabilities...")
    rejection_counter = Counter(0, 0, 0, 0, 0, 0)

    for i in 1:n_iterations

        if i % 100000 == 0
            println("$α_str\t $Δt_str\t  $i / $n_iterations")
            flush(stdout)
        end

        if scheme == "RMHMC"
            one_step_RMHMC!(q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, inv_sqrt_D, α, κ_α, Δt, rejection_counter)
        elseif scheme == "RMGHMC"
            one_step_RMGHMC!(q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, α, κ_α, Δt, γ, rejection_counter)
        end

    end

    proba_forward_no_convergence_momenta_update = rejection_counter.forward_no_convergence_momenta_update / n_iterations

    proba_forward_no_convergence_position_update = rejection_counter.forward_no_convergence_position_update / n_iterations

    proba_backward_no_convergence_momenta_update = rejection_counter.backward_no_convergence_momenta_update / n_iterations

    proba_backward_no_convergence_position_update = rejection_counter.backward_no_convergence_position_update / n_iterations

    proba_no_reversibility = rejection_counter.no_reversibility / n_iterations
    proba_MH_rejection = rejection_counter.MH_rejection / n_iterations

    proba_global_rejection = proba_forward_no_convergence_momenta_update + proba_forward_no_convergence_position_update + proba_backward_no_convergence_momenta_update + proba_backward_no_convergence_position_update + proba_no_reversibility + proba_MH_rejection

    println("
    proba_forward_no_convergence_momenta_update\t $proba_forward_no_convergence_momenta_update
    proba_forward_no_convergence_position_update\t $proba_forward_no_convergence_position_update
    proba_backward_no_convergence_momenta_update\t $proba_backward_no_convergence_momenta_update
    proba_backward_no_convergence_position_update\t $proba_backward_no_convergence_position_update
    proba_no_reversibility\t $proba_no_reversibility
    proba_MH_rejection\t $proba_MH_rejection
    proba_global_rejection\t $proba_global_rejection
    ")


    mkpath(dir_path)

    save_object(dir_path * "proba_forward_no_convergence_momenta_update.jld2", [proba_forward_no_convergence_momenta_update])
    save_object(dir_path * "proba_forward_no_convergence_position_update.jld2", [proba_forward_no_convergence_position_update])
    save_object(dir_path * "proba_backward_no_convergence_momenta_update.jld2", [proba_backward_no_convergence_momenta_update])
    save_object(dir_path * "proba_backward_no_convergence_position_update.jld2", [proba_backward_no_convergence_position_update])
    save_object(dir_path * "proba_no_reversibility.jld2", [proba_no_reversibility])
    save_object(dir_path * "proba_MH_rejection.jld2", [proba_MH_rejection])
    save_object(dir_path * "proba_global_rejection.jld2", [proba_global_rejection])

end


function compute_rejection_probabilities(scheme, Δt_range, n_iterations, α; γ=1.0)

    @threads for Δt in Δt_range
        rejection_probabilities_Δt(scheme, n_iterations, Δt, α; γ=γ)
    end
end



scheme = "RMHMC"
α = 0.6

if scheme == "RMHMC"
    Δt_min = 3e-2
    Δt_max = 1e-1
    small_Δt_values = [1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2, 2e-2]
    # α_values = [0.6,0.8,1.0,1.5]
elseif scheme == "RMGHMC"
    Δt_min = 5e-3
    Δt_max = 5e-2
    small_Δt_values = [1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3, 2e-3, 3e-3, 4e-3]
    # α_values = [0.6, 1.0, 1.4, 2.0]
end


Δt_values = 10.0 .^ LinRange(log(Δt_min) / log(10), log(Δt_max) / log(10), 16)
Δt_values = [small_Δt_values; Δt_values]

γ = 1.0
n_iterations = Int(1e7)

println("Running for $scheme with α=$α")
compute_rejection_probabilities(scheme, Δt_values, n_iterations, α; γ=γ)