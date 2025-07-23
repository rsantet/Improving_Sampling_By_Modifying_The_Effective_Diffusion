using CairoMakie

function plot_DW()
    RR = LinRange(0, box_length / 2, 1000)

    fig, ax, pl = lines(RR, DW.(RR), label="DW")
    ax.xlabel = L"r"
    ax.ylabel = L"V_{\mathrm{DW}}"
    ax.xlabelsize = 20
    ax.ylabelsize = 20
    ax.xticklabelsize = 15
    ax.yticklabelsize = 15
    # lines!(RR, ∇DW.(RR), label="∇DW")
    # axislegend()

    img_dir = "img/potentials/"
    mkpath(img_dir)
    save(img_dir * "DW.png", fig)


end

plot_DW()

function plot_WCA()

    RR = LinRange(r0 / 2, r0, 1000)

    fig, ax, pl = lines(RR, WCA.(RR), label="WCA")
    ax.xlabel = L"r"
    ax.ylabel = L"V_{\mathrm{WCA}}"
    ax.xlabelsize = 20
    ax.ylabelsize = 20
    ax.xticklabelsize = 15
    ax.yticklabelsize = 15
    # lines!(RR, ∇WCA.(RR), label="∇WCA")
    # axislegend()

    img_dir = "img/potentials/"
    mkpath(img_dir)
    save(img_dir * "WCA.png", fig)

end

plot_WCA()