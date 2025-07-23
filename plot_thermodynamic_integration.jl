using CairoMakie
using JLD2
using LaTeXStrings

include("input.jl")

function plot_TI()

    MF = load_object("data/TI/mean_force.jld2")
    FE = load_object("data/TI/free_energy.jld2")
    z_values = load_object("data/TI/z_values.jld2")
    fig, ax, pl = lines(z_values, MF, linewidth=3)
    pl.dpi = 300
    ax.xlabel = L"z"
    ax.ylabel = "Mean force"
    ax.xlabelsize = 35
    ax.ylabelsize = 30
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.limits = (-0.27, 1.23, -25, 25)

    # ax2 = Axis(fig, bbox=BBox(202,452,200,400))
    # ax2.xlabel=L"z"
    # ax2.ylabel="Mean force"
    # ax2.xlabelsize = 15
    # ax2.ylabelsize = 15
    # ax2.xticklabelsize = 10
    # ax2.yticklabelsize = 10
    # ax2.limits=(-0.27,1.23,-25,25)
    # ax2.xgridvisible=false
    # ax2.ygridvisible=false
    # lines!(ax2, z_values, MF)

    img_dir = "img/TI/"
    mkpath(img_dir)
    save(img_dir * "mean_force.png", fig)

    z_values_FE = [z_min + Δz * i for i in 0:nz]
    fig, ax, pl = lines(z_values_FE, FE, linewidth=3)
    pl.dpi = 300
    ax.xlabel = L"z"
    ax.ylabel = "Free energy"
    ax.xlabelsize = 35
    ax.ylabelsize = 30
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.limits = (-0.27, 1.23, -0.1, 3.5)

    # ax2 = Axis(fig, bbox=BBox(202,452,200,400))
    # ax2.xlabel=L"z"
    # ax2.ylabel="Free energy"
    # ax2.xlabelsize = 15
    # ax2.ylabelsize = 15
    # ax2.xticklabelsize = 10
    # ax2.yticklabelsize = 10
    # ax2.limits=(-0.27,1.23,-0.1,3.5)
    # ax2.xgridvisible=false
    # ax2.ygridvisible=false
    # lines!(ax2, z_values_FE, FE)

    save(img_dir * "free_energy.png", fig)

end

plot_TI()