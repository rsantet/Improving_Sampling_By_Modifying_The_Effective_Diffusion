include("input.jl")

using CairoMakie, Formatting

function plot_loggers(α, Δt, FE_plot_indices, MF_plot_indices)


    α_str = round(α, sigdigits=5)
    Δt_str = round(Δt, sigdigits=5)


    data_dir = "data/adaptive_MALA/$α_str/$Δt_str/"
    data_dir_fe = data_dir * "FE/"
    data_dir_mf = data_dir * "MF/"

    n_loggers = load_object(data_dir * "n_loggers.jld2")[1]
    size_loggers = load_object(data_dir * "size_loggers.jld2")[1]
    XX = collect(1:size_loggers) * n_loggers

    rc_values = load_object(data_dir * "rc_values.jld2")

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Number of iterations", ylabel=L"z", xlabelsize=25, ylabelsize=25, xticklabelsize=20, yticklabelsize=20)
    pltobj = lines!(XX, rc_values, label="")

    img_dir = "img/adaptive_MALA/$α_str/$Δt_str/"
    mkpath(img_dir)
    save(img_dir * "rc_values.png", fig)

    fig = Figure()
    ax = Axis(fig[1, 1], xlabelsize=25, ylabelsize=25, xticklabelsize=15, yticklabelsize=20, ylabel="Free energy", xlabel=L"z")
    pltobj = hist!(rc_values, normalization=:pdf, bins=80)
    # axislegend()
    save(img_dir * "rc_values_hist.png", fig)

    # Free Energy
    fig = Figure()
    ax = Axis(fig[1, 1], xlabelsize=25, ylabelsize=25, xticklabelsize=15, yticklabelsize=20, ylabel="Free Energy", xlabel=L"z")
    for i in FE_plot_indices
        FE_data = load_object(data_dir_fe * "$i.jld2")
        lines!(ax, z_values, FE_data, label=format(i, commas=true))
    end
    fig[1, 2] = Legend(fig, ax, "Iteration", titlesize=15, labelsize=15)
    save(img_dir * "FE.png", fig)

    # Mean force
    fig = Figure()
    ax = Axis(fig[1, 1], xlabelsize=25, ylabelsize=25, xticklabelsize=15, yticklabelsize=20, ylabel="Mean Force", xlabel=L"z")
    for i in MF_plot_indices
        MF_data = load_object(data_dir_mf * "$i.jld2")
        MF_data[MF_data[:, 1].==0, 2] .= 0.0
        data = MF_data[:, 2] ./ MF_data[:, 1]
        lines!(ax, z_values, data, label=format(i, commas=true))
    end
    fig[1, 2] = Legend(fig, ax, "Iteration", titlesize=15, labelsize=15)
    save(img_dir * "MF.png", fig)

end


α = 1.5
Δt = 0.002470885772418033
FE_plot_indices = [4000, 10_000, 20_000, 30_000, 40_000]
MF_plot_indices = [4000, 10_000, 20_000, 30_000, 40_000]

plot_loggers(α, Δt, FE_plot_indices, MF_plot_indices)