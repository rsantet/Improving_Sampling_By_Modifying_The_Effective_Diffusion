using JLD2
using CairoMakie
using LaTeXStrings
# using GLM
include("input.jl")

function plot_rejection_probabilities(scheme, Δt_range, n_iterations, α)

    α_str = round(α, sigdigits=5)

    MH_rejection_probabilities = []
    no_backward_convergence_probabilities_momenta = []
    no_backward_convergence_probabilities_position = []
    no_forward_convergence_probabilities_momenta = []
    no_forward_convergence_probabilities_position = []
    no_reversibility_probabilities = []
    global_rejection = []

    Δt_values_got = []


    for Δt in Δt_range
        Δt_str = round(Δt, sigdigits=5)
        dir_path = "data/$scheme/rejection_probabilities/$n_iterations/$α_str/$Δt_str/"

        if !isfile(dir_path * "proba_MH_rejection.jld2")
            continue
        end

        push!(Δt_values_got, Δt)

        proba_MH_rejection = load_object(dir_path * "proba_MH_rejection.jld2")
        proba_backward_no_convergence_momenta_update = load_object(dir_path * "proba_backward_no_convergence_momenta_update.jld2")
        proba_backward_no_convergence_position_update = load_object(dir_path * "proba_backward_no_convergence_position_update.jld2")
        proba_forward_no_convergence_momenta_update = load_object(dir_path * "proba_forward_no_convergence_momenta_update.jld2")
        proba_forward_no_convergence_position_update = load_object(dir_path * "proba_forward_no_convergence_position_update.jld2")
        proba_no_reversibility = load_object(dir_path * "proba_no_reversibility.jld2")
        proba_global_rejection = load_object(dir_path * "proba_global_rejection.jld2")

        push!(MH_rejection_probabilities, proba_MH_rejection[1])
        push!(no_backward_convergence_probabilities_momenta, proba_backward_no_convergence_momenta_update[1])
        push!(no_backward_convergence_probabilities_position, proba_backward_no_convergence_position_update[1])
        push!(no_forward_convergence_probabilities_momenta, proba_forward_no_convergence_momenta_update[1])
        push!(no_forward_convergence_probabilities_position, proba_forward_no_convergence_position_update[1])
        push!(no_reversibility_probabilities, proba_no_reversibility[1])
        push!(global_rejection, proba_global_rejection[1])

        println("$scheme\t $α_str \t $Δt_str")
        println("Forward Momenta:\t $(proba_forward_no_convergence_momenta_update[1])")
        println("Forward Position:\t $(proba_forward_no_convergence_position_update[1])")
        println("Backward Momenta:\t $(proba_backward_no_convergence_momenta_update[1])")
        println("Backward Position:\t $(proba_backward_no_convergence_position_update[1])")
        println("No Rev:\t $(proba_no_reversibility[1])")
        println("MH:\t $(proba_MH_rejection[1])")
        println("Global:\t $(proba_global_rejection[1])")
        println()

    end

    if length(Δt_values_got) == 0
        println("$scheme \t No values for α=$α")
        return
    end

    MH_rejection_probabilities_mask = MH_rejection_probabilities .> 0
    no_backward_convergence_probabilities_momenta_mask = no_backward_convergence_probabilities_momenta .> 0
    no_backward_convergence_probabilities_position_mask = no_backward_convergence_probabilities_position .> 0
    no_forward_convergence_probabilities_momenta_mask = no_forward_convergence_probabilities_momenta .> 0
    no_forward_convergence_probabilities_position_mask = no_forward_convergence_probabilities_position .> 0
    no_reversibility_probabilities_mask = no_reversibility_probabilities .> 0
    mask = global_rejection .> 0

    fig = Figure()
    ax = Axis(fig[1, 1], xlabelsize=20, ylabelsize=20, xticklabelsize=15, yticklabelsize=15, xlabel=L"\Delta t", ylabel="Rejection probability", xscale=log10, yscale=log10, xgridvisible=false, ygridvisible=false)

    if scheme == "RMGHMC" # α = 1.0
        bbox = BBox(85, 245, 260, 430)
    elseif scheme == "RMHMC" # α = 0.8
        bbox = BBox(85, 245, 280, 430)
    end

    ax2 = Axis(fig, bbox=bbox, xlabelsize=10, ylabelsize=10, xticklabelsize=7, yticklabelsize=7, xlabel=L"\Delta t", ylabel="Rejection probability", xscale=log10, yscale=log10, yaxisposition=:right, xgridvisible=false, ygridvisible=false)
    ylims!(ax, low=1e-5)
    ylims!(ax2, low=1e-5)

    if sum(MH_rejection_probabilities_mask) > 0
        scatterlines!(ax, Δt_values_got[MH_rejection_probabilities_mask], MH_rejection_probabilities[MH_rejection_probabilities_mask], label="Metropolis-Hastings", color="grey", marker=:diamond, markersize=10)
        scatterlines!(ax2, Δt_values_got[MH_rejection_probabilities_mask][end-15:end], MH_rejection_probabilities[MH_rejection_probabilities_mask][end-15:end], label="Metropolis-Hastings", color="grey", marker=:diamond, markersize=10)
    end
    if sum(no_backward_convergence_probabilities_momenta_mask) > 0
        scatterlines!(ax, Δt_values_got[no_backward_convergence_probabilities_momenta_mask], no_backward_convergence_probabilities_momenta[no_backward_convergence_probabilities_momenta_mask], label="Backward (momenta)", color="lightgreen", marker=:utriangle, markersize=10)
        scatterlines!(ax2, Δt_values_got[no_backward_convergence_probabilities_momenta_mask], no_backward_convergence_probabilities_momenta[no_backward_convergence_probabilities_momenta_mask], label="Backward (momenta)", color="lightgreen", marker=:utriangle, markersize=10)
    end
    if sum(no_backward_convergence_probabilities_position) > 0
        scatterlines!(ax, Δt_values_got[no_backward_convergence_probabilities_position_mask], no_backward_convergence_probabilities_position[no_backward_convergence_probabilities_position_mask], label="Backward (position)", color="darkgreen", marker=:dtriangle, markersize=10)
        if scheme == "RMGHMC" # α = 1.0
            scatterlines!(ax2, Δt_values_got[no_backward_convergence_probabilities_position_mask], no_backward_convergence_probabilities_position[no_backward_convergence_probabilities_position_mask], label="Backward (position)", color="darkgreen", marker=:dtriangle, markersize=10)
        elseif scheme == "RMHMC" #α =0.8
            scatterlines!(ax2, Δt_values_got[no_backward_convergence_probabilities_position_mask][2:end], no_backward_convergence_probabilities_position[no_backward_convergence_probabilities_position_mask][2:end], label="Backward (position)", color="darkgreen", marker=:dtriangle, markersize=10)
        end
    end
    if sum(no_forward_convergence_probabilities_momenta_mask) > 0
        scatterlines!(ax, Δt_values_got[no_forward_convergence_probabilities_momenta_mask], no_forward_convergence_probabilities_momenta[no_forward_convergence_probabilities_momenta_mask], label="Forward (momenta)", color="lightblue", marker=:rtriangle, markersize=10)
        if scheme == "RMGHMC" # α=1.0
            scatterlines!(ax2, Δt_values_got[no_forward_convergence_probabilities_momenta_mask][4:end], no_forward_convergence_probabilities_momenta[no_forward_convergence_probabilities_momenta_mask][4:end], label="Forward (momenta)", color="lightblue", marker=:rtriangle, markersize=10)
        elseif scheme == "RMHMC" # α=0.8
            scatterlines!(ax2, Δt_values_got[no_forward_convergence_probabilities_momenta_mask][9:end], no_forward_convergence_probabilities_momenta[no_forward_convergence_probabilities_momenta_mask][9:end], label="Forward (momenta)", color="lightblue", marker=:rtriangle, markersize=10)
        end
    end
    if sum(no_forward_convergence_probabilities_position_mask) > 0
        scatterlines!(ax, Δt_values_got[no_forward_convergence_probabilities_position_mask], no_forward_convergence_probabilities_position[no_forward_convergence_probabilities_position_mask], label="Forward (position)", color="darkblue", marker=:ltriangle, markersize=10)
        if scheme == "RMGHMC" # α=1.0
            scatterlines!(ax2, Δt_values_got[no_forward_convergence_probabilities_position_mask][2:end], no_forward_convergence_probabilities_position[no_forward_convergence_probabilities_position_mask][2:end], label="Forward (position)", color="darkblue", marker=:ltriangle, markersize=10)
        elseif scheme == "RMHMC" # α=0.8
            scatterlines!(ax2, Δt_values_got[no_forward_convergence_probabilities_position_mask][6:end], no_forward_convergence_probabilities_position[no_forward_convergence_probabilities_position_mask][6:end], label="Forward (position)", color="darkblue", marker=:ltriangle, markersize=10)
        end
    end
    if sum(no_reversibility_probabilities_mask) > 0
        scatterlines!(ax, Δt_values_got[no_reversibility_probabilities_mask], no_reversibility_probabilities[no_reversibility_probabilities_mask], label="Reversibility", color="red", marker=:pentagon, markersize=10)
        if scheme == "RMGHMC" # α = 1.0
            scatterlines!(ax2, Δt_values_got[no_reversibility_probabilities_mask], no_reversibility_probabilities[no_reversibility_probabilities_mask], label="Reversibility", color="red", marker=:pentagon, markersize=10)
        elseif scheme == "RMHMC" # α = 0.8
            scatterlines!(ax2, Δt_values_got[no_reversibility_probabilities_mask][4:end], no_reversibility_probabilities[no_reversibility_probabilities_mask][4:end], label="Reversibility", color="red", marker=:pentagon, markersize=10)
        end
    end
    if sum(mask) > 0
        scatterlines!(ax, Δt_values_got[mask], global_rejection[mask], label="Global", color="black", marker=:star4, markersize=10)
        scatterlines!(ax2, Δt_values_got[mask][end-15:end], global_rejection[mask][end-15:end], label="Global", color="black", marker=:star4, markersize=10)
        # lines!(ax, Δt_values_got[MH_rejection_probabilities_mask], Δt_values_got[MH_rejection_probabilities_mask] .^ 3 .* 6000, label=L"\Delta t^3")
        # lines!(ax2, Δt_values_got[MH_rejection_probabilities_mask][end-15:end], Δt_values_got[MH_rejection_probabilities_mask][end-15:end] .^ 3 .* 6000, label=L"\Delta t^3")
        # lines!(ax, Δt_values_got[MH_rejection_probabilities_mask], Δt_values_got[MH_rejection_probabilities_mask] .^ (3 / 2) .* 6000, label=L"\Delta t^{3/2}")
    end



    fig[1, 2] = Legend(fig, ax, titlesize=20, labelsize=15)
    img_dir = "img/$scheme/rejection_probabilities/$n_iterations/$α_str/"
    mkpath(img_dir)

    save(img_dir * "rejection_probabilities.png", fig)

end


γ = 1.0
n_iterations = Int(1e7)
scheme = "RMGHMC"

if scheme == "RMHMC"
    Δt_min = 3e-2
    Δt_max = 1e-1
    small_Δt_values = [1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2, 2e-2]
    α_values = [0.6, 0.8, 1.0, 1.5]
elseif scheme == "RMGHMC"
    Δt_min = 5e-3
    Δt_max = 5e-2
    small_Δt_values = [1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3, 2e-3, 3e-3, 4e-3]
    α_values = [0.6, 1.0, 1.4, 2.0]
end


Δt_values = 10.0 .^ LinRange(log(Δt_min) / log(10), log(Δt_max) / log(10), 16)
Δt_values = [small_Δt_values; Δt_values]

for α in α_values
    println("Running for $scheme with α=$α")
    plot_rejection_probabilities(scheme, Δt_values, n_iterations, α)
end