include("input.jl")
include("schemes.jl")

function run_adaptive_MALA(initial_z, α, n_iterations, Δt; save::Bool=false, n_loggers=10_000)

    α_str = round(α, sigdigits=5)
    Δt_str = round(Δt, sigdigits=5)

    q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, inv_D, div_D, MeanForce, FreeEnergy, MeanForceUpdate, κ_α, sqrt_κ_α = initialize_adaptive_OL(initial_z, α)

    rejection_counter = Counter(0, 0, 0, 0, 0, 0)

    if save
        if !(typeof(n_loggers) <: Int)
            println("Use an Int for n_loggers when saving")
            return
        end
        data_dir = "data/adaptive_MALA/$α_str/$Δt_str/"
        mkpath(data_dir)
        data_dir_fe = data_dir * "FE/"
        data_dir_mf = data_dir * "MF/"
        mkpath(data_dir_fe)
        mkpath(data_dir_mf)
        size_loggers = floor(Int, n_iterations / n_loggers)
        save_object(data_dir * "n_loggers.jld2", [n_loggers])
        save_object(data_dir * "size_loggers.jld2", [size_loggers])

        counter_loggers = 0
        traj = Vector{Vector{SVector{2,Float64}}}(undef, size_loggers)
        rc_values = Vector{Float64}(undef, size_loggers)
    end

    for i in 1:n_iterations

        if i % n_mean_force_update == 0
            update_mean_force!(MeanForceUpdate, MeanForce)
            update_free_energy!(MeanForce, FreeEnergy)
            κ_α, sqrt_κ_α = update_κ_α(FreeEnergy, α)
        end

        if i % 10_000 == 0
            z = ξ(q, unit_vector, distance)
            println("α = $α_str\t Δt=$Δt_str\t $i / $n_iterations\t CV VALUE\t $z")
        end

        z_value = one_step_adaptive_MALA!(q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, div_D, inv_D, Δt, α, κ_α, sqrt_κ_α, MeanForceUpdate, MeanForce, FreeEnergy, rejection_counter)

        if save
            if i % n_loggers == 0
                counter_loggers += 1
                traj[counter_loggers] = deepcopy(q)
                rc_values[counter_loggers] = z_value

                save_object(data_dir_fe * "$i.jld2", FreeEnergy)
                save_object(data_dir_mf * "$i.jld2", MeanForce)

                # println("κ_α is $κ_α")
                # println(FreeEnergy)
            end
        end

    end

    proba_MH_rejection = rejection_counter.MH_rejection / n_iterations


    println("
    proba_MH_rejection\t $proba_MH_rejection
    ")


    if save
        save_object(data_dir * "traj.jld2", traj)
        save_object(data_dir * "rc_values.jld2", rc_values)

        save_object(data_dir * "proba_MH_rejection.jld2", [proba_MH_rejection])
    end

end

initial_z = 0.0
α = 1.5
n_iterations = 40_000
Δt = 0.002470885772418033
n_loggers = 500
run_adaptive_MALA(initial_z, α, n_iterations, Δt; save=true, n_loggers=n_loggers)