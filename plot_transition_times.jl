include("input.jl")
using CairoMakie
using LaTeXStrings
using Statistics


function plot_transition_times(scheme, Δt_values, α_values)

    println("Plotting for $scheme")

    ##################### Preparing Figures #####################

    fig_tt = Figure()
    ax_tt = Axis(fig_tt[1, 1], xlabelsize=20, ylabelsize=20, xticklabelsize=15, yticklabelsize=15, xlabel=L"\Delta t", xscale=log10, yscale=log10)

    fig_n_it = Figure()
    ax_n_it = Axis(fig_n_it[1, 1], xlabelsize=25, ylabelsize=25, xticklabelsize=20, yticklabelsize=20, xlabel=L"\Delta t", ylabel="Mean number of iterations", xscale=log10, yscale=log10)

    fig_min_vals = Figure()
    ax_min_vals = Axis(fig_min_vals[1, 1], xlabelsize=25, ylabelsize=25, xticklabelsize=20, yticklabelsize=20, xlabel=L"\alpha", ylabel="Mean number of iterations")


    ##################### Nonconstant diffusions #####################

    min_vals = Float64[]
    min_vals_std = Float64[]
    min_vals_n_transitions = Float64[]
    argmin_vals = Float64[]
    α_min_vals = Float64[]


    if length(α_values) > 1
        palette = cgrad(:rainbow, LinRange(0, 1, length(α_values)))
    else
        palette = [:red]
    end

    colormap_matrix = zeros(length(α_values), length(Δt_values))

    for (idx_α, α) in enumerate(α_values)

        tt_vec_means = Float64[]
        tt_vec_stds = Float64[]
        tt_n_transitions = Float64[]
        n_iterations_vec_means = Float64[]
        n_iterations_vec_stds = Float64[]
        n_iterations_n_transitions = Float64[]
        Δt_values_got = Float64[]
        min_val = Inf
        min_val_std = Inf
        min_val_n_transitions = Inf
        argmin_val = 0.0

        for (idx_Δt, Δt) in enumerate(Δt_values)
            data_dir = "data/$scheme/transition_times/$(round(α, sigdigits=5))/$(round(Δt, sigdigits=5))/"

            tt_fp = data_dir * "tt.jld2"
            n_iterations_fp = data_dir * "n_iterations.jld2"

            if isfile(tt_fp) && isfile(n_iterations_fp)

                tt = load_object(tt_fp)
                n_iterations = load_object(n_iterations_fp)
                len_tt = length(tt)
                len_n_it = length(n_iterations)
                # println(Δt)
                # println(α)
                # println(length(tt))
                # @assert length(tt) == 10000



                tt_mean = mean(tt)
                tt_std = std(tt)
                n_it_mean = mean(n_iterations)
                n_it_std = std(n_iterations)


                push!(tt_vec_means, tt_mean)
                push!(tt_vec_stds, tt_std)
                push!(n_iterations_vec_means, n_it_mean)
                push!(n_iterations_vec_stds, n_it_std)
                push!(tt_n_transitions, len_tt)
                push!(n_iterations_n_transitions, len_n_it)

                colormap_matrix[idx_α, idx_Δt] = n_it_mean

                push!(Δt_values_got, Δt)

                if n_it_mean < min_val
                    min_val = n_it_mean
                    min_val_std = n_it_std
                    min_val_n_transitions = len_n_it
                    argmin_val = Δt
                end

                # println("α = $(round(α, sigdigits=5))\t Δt = $(round(Δt, sigdigits=5))\t Found $(len_n_it)/$len_tt transitions")

            else

                println("α = $(round(α, sigdigits=5))\t Missing Δt = $(round(Δt, sigdigits=5))")

            end

        end

        if min_val < Inf
            push!(min_vals, min_val)
            push!(min_vals_std, min_val_std)
            push!(min_vals_n_transitions, min_val_n_transitions)
            push!(α_min_vals, α)
            push!(argmin_vals, argmin_val)
        end

        ##################### PLOTTING NONCONSTANT DIFFUSIONS #####################

        if length(tt_vec_means) > 0 && length(n_iterations_vec_means) > 0

            tt_CI_error = 1.96 .* tt_vec_stds ./ sqrt.(tt_n_transitions)

            n_iterations_CI_error = 1.96 .* n_iterations_vec_stds ./ sqrt.(n_iterations_n_transitions)

            # errorbars!(ax_tt, Δt_values_got, tt_vec_means, tt_CI_error, whiskerwidth=10, color=palette[idx_α])
            # errorbars!(ax_n_it, Δt_values_got, n_iterations_vec_means, n_iterations_CI_error, whiskerwidth=10, color=palette[idx_α])

            scatter!(ax_n_it, Δt_values_got, n_iterations_vec_means, label=L"\alpha=%$(round(α, digits=5))", color=palette[idx_α])
            # scatter!(ax_n_it, Δt_values_got, n_iterations_vec_means, label="CV Diffusion", color=palette[idx_α], markersize=10)
            scatter!(ax_tt, Δt_values_got, tt_vec_means, label=L"\alpha=%$(round(α, digits=5))", color=palette[idx_α])

        else

            # println("No data for α = $(round(α, sigdigits=5))")

        end
    end

    ##################### Constant diffusion #####################

    constant_tt_vec_means = Float64[]
    constant_tt_vec_stds = Float64[]
    constant_tt_n_transitions = Float64[]
    constant_n_iterations_vec_means = Float64[]
    constant_n_iterations_vec_stds = Float64[]
    constant_n_iterations_n_transitions = Float64[]
    constant_Δt_values_got = Float64[]
    constant_min_val = Inf
    constant_min_val_std = Inf
    constant_min_val_n_transitions = Inf
    constant_argmin_val = 0.0

    for (idx_Δt, Δt) in enumerate(Δt_values)
        data_dir = "data/$scheme/transition_times/constant_diffusion/$(round(Δt, sigdigits=5))/"

        tt_fp = data_dir * "tt.jld2"
        n_iterations_fp = data_dir * "n_iterations.jld2"

        if isfile(tt_fp) && isfile(n_iterations_fp)

            tt = load_object(tt_fp)
            n_iterations = load_object(n_iterations_fp)
            len_tt = length(tt)
            len_n_it = length(n_iterations)
            # println(Δt)
            # println(α)
            # println(length(tt))
            # @assert length(tt) == 10000

            tt_mean = mean(tt)
            tt_std = std(tt)
            n_it_mean = mean(n_iterations)
            n_it_std = std(n_iterations)



            push!(constant_tt_vec_means, tt_mean)
            push!(constant_tt_vec_stds, tt_std)
            push!(constant_n_iterations_vec_means, n_it_mean)
            push!(constant_n_iterations_vec_stds, n_it_std)
            push!(constant_tt_n_transitions, len_tt)
            push!(constant_n_iterations_n_transitions, len_n_it)

            push!(constant_Δt_values_got, Δt)

            if n_it_mean < constant_min_val
                constant_min_val = n_it_mean
                constant_min_val_std = n_it_std
                constant_min_val_n_transitions = len_n_it
                constant_argmin_val = Δt
            end

            # println("Constant diffusion\t Δt = $(round(Δt, sigdigits=5))\t Found $(len_n_it)/$len_tt transitions")
        else
            println("Constant diffusion\t Missing Δt = $(round(Δt, sigdigits=5))")
        end
    end


    ##################### PLOTTING CONSTANT DIFFUSION #####################

    if length(constant_tt_vec_means) > 0 && length(constant_n_iterations_vec_means) > 0

        constant_tt_CI_error = 1.96 .* constant_tt_vec_stds ./ sqrt.(constant_tt_n_transitions)

        constant_n_iterations_CI_error = 1.96 .* constant_n_iterations_vec_stds ./ sqrt.(constant_n_iterations_n_transitions)

        # errorbars!(ax_tt, constant_Δt_values_got, constant_tt_vec_means, constant_tt_CI_error, whiskerwidth=10, color="black")
        # errorbars!(ax_n_it, constant_Δt_values_got, constant_n_iterations_vec_means, constant_n_iterations_CI_error, whiskerwidth=10, color="black")

        scatter!(ax_n_it, constant_Δt_values_got, constant_n_iterations_vec_means, label="Constant diffusion", color="black", markersize=10)
        scatter!(ax_tt, constant_Δt_values_got, constant_tt_vec_means, label="Constant diffusion", color="black")

    else
        # println("No data for Constant diffusion")
    end

    ##################### LEGENDS #####################

    img_dir = "img/$scheme/transition_times/"
    mkpath(img_dir)

    # Plot behavior in the limit \Delta t\to 0
    # if scheme == "RMGHMC"
    #     lines!(ax_n_it, Δt_values[1:6], Δt_values[1:6] .^ (-1) .* 8, label=L"Δt^{-1}")
    # elseif scheme == "RMHMC"
    #     lines!(ax_n_it, Δt_values[1:6], Δt_values[1:6] .^ (-2) .* 10, label=L"Δt^{-2}")
    # elseif scheme == "MALA" || scheme == "adaptive_MALA"
    #     lines!(ax_n_it, Δt_values[1:6], Δt_values[1:6] .^ (-1) .*8, label=L"Δt^{-1}")
    # end

    # fig_n_it[1, 2] = Legend(fig_n_it, ax_n_it, labelsize=15)#, nbanks=3)
    fig_tt[1, 2] = Legend(fig_tt, ax_tt, labelsize=15)#, nbanks=3)
    # lines!(ax_n_it,Δt_values[1:6], Δt_values[1:6] .^(-2), label=L"Δt^{-2}")
    # lines!(ax_tt,Δt_values[1:6], Δt_values[1:6] .^(-1), label=L"Δt^{-2}")

    axislegend(ax_tt, position=:lt)
    if scheme == "MALA" || scheme == "adaptive_MALA"
        axislegend(ax_n_it, position=:lt, labelsize=20)
    else
        axislegend(ax_n_it, position=:rt, labelsize=20)
    end
    save(img_dir * "transition_times.png", fig_tt)
    save(img_dir * "n_iterations.png", fig_n_it)


    ##################### MINIMUM VALUES OBTAINED #####################

    fig_colormap = Figure()
    ax_colormap = Axis(fig_colormap[1, 1], xlabelsize=20, ylabelsize=20, xlabel=L"\alpha", ylabel=L"\Delta t", yscale=log10)
    hm = heatmap!(ax_colormap, α_values, Δt_values, colormap_matrix, colormap=:viridis)
    Colorbar(fig_colormap[1, 2], hm, label="Mean number of iterations")
    save(img_dir * "colormap.png", fig_colormap)

    if length(min_vals) > 0 && length(α_min_vals) > 0

        # scatter!(ax_min_vals, α_min_vals, min_vals, color=argmin_vals, colormap=Makie.Categorical(:viridis))
        # unique_vals = unique(argmin_vals)
        # colorrange = [minimum(unique_vals), maximum(unique_vals)]
        # cbar = Colorbar(fig_min_vals[1, 2], colorrange=colorrange, colormap=cgrad(:viridis, length(unique_vals), categorical=true), tickformat="{:.2e}", ticks=unique_vals)
        # println(typeof(cbar.ticks))
        # fig_min_vals[1, 2] = Legend(fig_min_vals, ax_min_vals, framevisible=false, labelsize=12)
        scatter!(ax_min_vals, α_min_vals, min_vals, label=nothing)
        errorbars!(ax_min_vals, α_min_vals, min_vals, 1.96 .* min_vals_std ./ sqrt.(min_vals_n_transitions), whiskerwidth=10, color=:black, label=nothing)

        # Constant diffusion
        constant_min_val_CI_error = 1.96 .* constant_min_val_std ./ sqrt.(constant_min_val_n_transitions)
        hlines!(ax_min_vals, constant_min_val, 0, 1, label="Constant diffusion", color=:red, linestyle=:dash)
        band_lower = (constant_min_val - constant_min_val_CI_error) * ones(length(α_min_vals))
        band_upper = (constant_min_val + constant_min_val_CI_error) * ones(length(α_min_vals))
        band!(ax_min_vals, α_min_vals, band_lower, band_upper, color=("red", 0.1))
        if scheme == "MALA" || scheme == "adaptive_MALA"
            axislegend(ax_min_vals, position=:rt, labelsize=20)
        else
            axislegend(ax_min_vals, position=:lt, labelsize=20)
        end
        save(img_dir * "min_vals.png", fig_min_vals)

        for i in eachindex(min_vals)
            println("α=$(α_values[i])\t Optimal time step $(argmin_vals[i])")
        end
        println(argmin_vals)
        println(min_vals)

        println("Value for Constant diffusion is $constant_min_val")
        argmin_val = argmin(min_vals)
        println("Value for α = $(round(α_min_vals[argmin_val], sigdigits=5)) is $(round(min_vals[argmin_val], sigdigits=5))")
    end

    println()
end



# schemes = ["adaptive_MALA", "MALA", "RMHMC", "RMGHMC"]
schemes = ["RMGHMC"]
α_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4]
γ = 1.0

for scheme in schemes
    if scheme == "RMGHMC"
        Δt_min = 5e-3
        Δt_max = 5e-2
        # α_values = [0.6, 1.0] # For n_iterations plots
    elseif scheme == "RMHMC"
        Δt_min = 3e-2
        Δt_max = 1e-1
        # α_values = [0.6, 0.8] # For n_iterations plots
    elseif scheme == "MALA"
        Δt_min = 5e-4
        Δt_max = 1e-2
        # α_values = [0.6, 1.4] # For n_iterations plots
    elseif scheme == "adaptive_MALA"
        Δt_min = 5e-4
        Δt_max = 1e-2
        # α_values = [0.6, 1.5] # For n_iterations plots
    end
    
    Δt_values = 10.0 .^ LinRange(log(Δt_min) / log(10), log(Δt_max) / log(10), 16)

    plot_transition_times(scheme, Δt_values, α_values)
end