include("input_TI.jl")

using Statistics, JLD2, LinearAlgebra

function constrained_overdamped_langevin!(q, F, unit_vector, distance, Δt, current_z)

    # Force computation
    compute_forces!(q, F, unit_vector, distance)

    # Unconstrained move
    for i in eachindex(q)
        q[i] += Δt * F[i] + sqrt(2 * Δt / β) * SVector(randn(), randn())
    end

    # Computing the Lagrange multiplier
    z = ξ(q, unit_vector, distance)
    lagrange = (current_z - z) * 2 * w^2
    factor = lagrange / w / (r1 + 2 * w * current_z)
    if abs(factor - 1) > TI_tolerance
        factor = 1 / (1 - factor)
    else
        println("LAGRANGE MULTIPLIER UNDEFINED")
        factor = 0.0
    end

    # Deal with periodic conditions first
    nb = round.((q[1] - q[2]) / box_length)
    q[1] -= nb * box_length

    # Compute new_positions
    old_coords_1 = copy(q[1])
    old_coords_2 = copy(q[2])
    q[1] = 0.5 * (1 + factor) * old_coords_1 + 0.5 * (1 - factor) * old_coords_2
    q[2] = 0.5 * (1 - factor) * old_coords_1 + 0.5 * (1 + factor) * old_coords_2

    # Get back to the original frame
    q[1] += nb * box_length

    # Compute the local mean force
    lmf = local_mean_force(q, F, unit_vector, distance)

    return lmf
end

function simulate!(n_iterations::Int, q, F, unit_vector, distance, Δt, current_z; save_lmf_values::Bool=true)

    if save_lmf_values
        lmf_values = Vector{Float64}(undef, n_iterations)
    end

    @inbounds for i in 1:n_iterations
        if i % 1_000_000 == 0
            println("$i / $n_iterations")
        end
        lmf_value = constrained_overdamped_langevin!(q, F, unit_vector, distance, Δt, current_z)
        if save_lmf_values
            lmf_values[i] = lmf_value
        end
    end

    if save_lmf_values
        return lmf_values
    end
end

function TI(Δt, n_thermalization_iterations::Int, n_simulation_iterations::Int, z_min, z_max, nz::Int)

    # Initialize the system with CV value at z_min
    q, F, unit_vector, distance = initialize_OL(z_min)

    z_values = [z_min + (z_max - z_min) / nz * (i - 1 / 2) for i in 1:nz]

    data_dir = "data/TI/"
    mkpath(data_dir)

    lmf_dir = data_dir * "lmf/"
    if isdir(lmf_dir) # Remove old values if needed
        rm(lmf_dir, recursive=true)
    end
    mkpath(lmf_dir)

    for current_z in z_values

        # Initial thermalization
        println("Thermalization for z=\t$current_z")
        simulate!(n_thermalization_iterations, q, F, unit_vector, distance, Δt, current_z; save_lmf_values=false)

        println("Simulation for z=\t$current_z")

        lmf_values = simulate!(n_simulation_iterations, q, F, unit_vector, distance, Δt, current_z; save_lmf_values=true)

        println("Mean Force for z=$current_z:\t ", mean(lmf_values))

        # Save
        lmf_fp = lmf_dir * "$current_z.jld2"
        save_object(lmf_fp, lmf_values)
    end

end

function construct_MF_FE(z_min, z_max, nz)
    data_dir = "data/TI/"
    data_lmf_dir = data_dir * "lmf/"
    z_values = [z_min + (z_max - z_min) / nz * (i - 1 / 2) for i in 1:nz]
    Δz = (z_max-z_min)/nz

    MF = zeros(nz)
    FE = zeros(nz+1)

    for (i,z) in enumerate(z_values)
        fp = data_lmf_dir * "$z.jld2"
        lmf_values = load_object(fp)
        MF[i] = mean(lmf_values)
    end

    for i in 1:nz
        FE[i+1] = FE[i]+ MF[i]*Δz
    end
    FE .-= minimum(FE)

    save_object(data_dir * "z_values.jld2", z_values)
    save_object(data_dir * "free_energy.jld2", FE)
    save_object(data_dir * "mean_force.jld2", MF)
end

z_min = -0.5
z_max = 1.5
nz = 100

Δt = 2.5e-5
n_thermalization_iterations = floor(Int, 1 / Δt)
n_simulation_iterations = floor(Int, 125 / Δt)

TI(Δt, n_thermalization_iterations, n_simulation_iterations, z_min ,z_max, nz)

construct_MF_FE(z_min,z_max,nz)