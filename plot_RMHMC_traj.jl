using CairoMakie
using JLD2


α = 0.
Δt = 1e-3

α_str = round(α,sigdigits=5)
Δt_str = round(Δt,sigdigits=5)
dir_path = "data/RMHMC/$α_str/$Δt_str/"
traj = load_object(dir_path * "traj.jld2")
z_values = load_object(dir_path * "z_values.jld2")
lines(1:n_iterations, z_values)
for i in eachindex(traj)
    periodize!(traj[i])
end
function plot_traj(traj, z_values)
    x, y = [traj[i][j][1] for i in eachindex(traj), j in eachindex(traj[1])], [traj[i][j][2] for i in eachindex(traj), j in eachindex(traj[1])]

    fig, ax, pl = scatter(x[1,3:end],y[1,3:end])
    pl2 = scatter!(x[1,1:2], y[1,1:2], color="red", markersize=10)
    ax.limits = (0,box_length, 0, box_length)
    pl3 = annotations!(["$(round(z_values[1],sigdigits=1))"],[traj[1][1]])

    record(fig, "data/RMHMC_trajectory.mp4", 2:100:100_000; framerate=30) do i
        pl[1] = x[i,3:end]
        pl[2] = y[i,3:end]
        pl2[1] = x[i,1:2]
        pl2[2] = y[i,1:2]
        pl3[1] = ["$(round(z_values[i],sigdigits=2))"]
        pl3[2] = [traj[i][1]]
    end
end

plot_traj(traj, z_values)