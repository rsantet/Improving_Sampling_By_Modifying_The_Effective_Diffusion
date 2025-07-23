include("input.jl")
include("schemes.jl")

function run_MALA(initial_z, α, n_iterations, Δt; save::Bool=false)

    α_str = round(α, sigdigits=5)
    Δt_str = round(Δt,sigdigits=5)

    q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, inv_D, div_D, κ_α = initialize_OL(initial_z, α)

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

        z_value = one_step_MALA!(q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, div_D, inv_D, Δt, α, κ_α, rejection_counter)

        if save
            push!(traj, deepcopy(q))
            push!(z_values, z_value)
        end

    end

    proba_MH_rejection = rejection_counter.MH_rejection / n_iterations


    println("
    proba_MH_rejection\t $proba_MH_rejection
    ")


    dir_path = "data/MALA/$α_str/$Δt_str/"
    mkpath(dir_path)

    save_object(dir_path * "proba_MH_rejection.jld2",[proba_MH_rejection])

    if save
        save_object(dir_path * "traj.jld2", traj)
        save_object(dir_path * "z_values.jld2", z_values)
    end

end

initial_z = 0.
α = 0.
n_iterations = Int(1e5)
Δt = 1e-4
run_MALA(initial_z, α, n_iterations, Δt; save=true)