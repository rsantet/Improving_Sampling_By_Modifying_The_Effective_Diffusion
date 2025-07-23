using CairoMakie

function plot_CV()

    function ξ_r(r)
        return (r - r1) / 2 / w
    end

    function ∇ξ_r(r)
        return 1 / 2 / w
    end

    RR = LinRange(0, box_length / 2, 1000)
    fig, ax, pl = lines(RR, ξ_r.(RR), label="DW")
    ax.xlabel = L"r"
    ax.ylabel = L"\xi"
    ax.xlabelsize = 20
    ax.ylabelsize = 20
    ax.xticklabelsize = 15
    ax.yticklabelsize = 15
    # lines!(RR, ∇ξ_r.(RR), label="∇DW")
    # axislegend()

    img_dir = "img/CV/"
    mkpath(img_dir)
    save(img_dir * "CV.png", fig)

end
plot_CV()