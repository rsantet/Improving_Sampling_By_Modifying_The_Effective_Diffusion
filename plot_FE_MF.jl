include("input.jl")

using CairoMakie

function plot_FE_MF_functions()

    ZZ = LinRange(z_min-0.25,z_max+0.25,1000)

    fig, ax, pl = lines(ZZ, mean_force.(ZZ))
    ax.xlabel = L"z"
    ax.ylabel = "Mean force"
    ax.xlabelsize = 20
    ax.ylabelsize = 20
    ax.xticklabelsize = 15
    ax.yticklabelsize = 15
    
    img_dir = "img/FE_MF/"
    mkpath(img_dir)
    save(img_dir * "mean_force.png", fig)

    fig, ax, pl = lines(ZZ, free_energy.(ZZ))
    ax.xlabel = L"z"
    ax.ylabel = "Free energy"
    ax.xlabelsize = 20
    ax.ylabelsize = 20
    ax.xticklabelsize = 15
    ax.yticklabelsize = 15
    
    img_dir = "img/FE_MF/"
    mkpath(img_dir)
    save(img_dir * "free_energy.png", fig)

end

plot_FE_MF_functions()