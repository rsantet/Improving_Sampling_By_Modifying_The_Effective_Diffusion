include("input.jl")
include("schemes.jl")

function run_RMHMC(initial_z, α, n_iterations, Δt; save::Bool=false)

    α_str = round(α, sigdigits=5)
    Δt_str = round(Δt,sigdigits=5)

    q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, inv_sqrt_D, κ_α = initialize_HMC(initial_z, α)

    rejection_counter = Counter(0, 0, 0, 0, 0, 0)

    if save
        traj = []
        z_values = []
    end

    for i in 1:n_iterations

        if i % 10000 == 0
            z = ξ(q, unit_vector, distance)
            println("α = $α_str\t Δt=$Δt_str\t $i / $n_iterations\t CV VALUE\t $z")
        end

        one_step_RMHMC!(q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, inv_sqrt_D, α, κ_α, Δt, rejection_counter)

        if save
            push!(traj, deepcopy(q))
            z = ξ(q, unit_vector, distance)
            push!(z_values, z)
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


    dir_path = "data/RMHMC/$α_str/$Δt_str/"
    mkpath(dir_path)

    save_object(dir_path * "proba_forward_no_convergence_momenta_update.jld2", [proba_forward_no_convergence_momenta_update])
    save_object(dir_path * "proba_forward_no_convergence_position_update.jld2", [proba_forward_no_convergence_position_update])
    save_object(dir_path * "proba_backward_no_convergence_momenta_update.jld2", [proba_backward_no_convergence_momenta_update])
    save_object(dir_path * "proba_backward_no_convergence_position_update.jld2", [proba_backward_no_convergence_position_update])
    save_object(dir_path * "proba_no_reversibility.jld2", [proba_no_reversibility])
    save_object(dir_path * "proba_MH_rejection.jld2",[proba_MH_rejection])
    save_object(dir_path * "proba_global_rejection.jld2",[proba_global_rejection])

    if save
        save_object(dir_path * "traj.jld2", traj)
        save_object(dir_path * "z_values.jld2", z_values)
    end

end

initial_z = 0.
α = 0.
n_iterations = Int(1e5)
Δt = 1e-2
run_RMHMC(initial_z, α, n_iterations, Δt; save=true)