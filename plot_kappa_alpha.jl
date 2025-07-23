include("input.jl")

using CairoMakie
using LaTeXStrings

# Free Energy barrier is 3.119875510330223

function plot_κ_α()

    ZZ = LinRange(0,2.5,1000)

    fig, ax, pl = lines(ZZ, κ.(ZZ))
    vlines!([1+log(effective_diffusion_squared)/3.119875510330223], 0,1, label=L"1+\frac{ln(\sigma^2)}{\beta h}", color="red")
    ax.xlabel = L"\alpha"
    ax.ylabel = L"\kappa_\alpha"
    ax.xlabelsize = 20
    ax.ylabelsize = 20
    ax.xticklabelsize = 15
    ax.yticklabelsize = 15
    
    img_dir = "img/"
    axislegend()
    mkpath(img_dir)
    save(img_dir * "kappa.png", fig)

end

plot_κ_α()