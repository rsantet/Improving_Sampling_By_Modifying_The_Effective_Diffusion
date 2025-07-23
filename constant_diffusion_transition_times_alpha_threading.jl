include("constant_diffusion_transition_times.jl")

using Base.Threads

n_transitions = 10_000

effective_diffusion_squared = 1.0

scheme = "MALA"

if scheme == "RMGHMC"
    Δt_min = 5e-3
    Δt_max = 5e-2
elseif scheme == "RMHMC"
    Δt_min = 3e-2
    Δt_max = 1e-1
elseif scheme == "MALA" || scheme == "adaptive_MALA"
    Δt_min = 5e-4
    Δt_max = 1e-2
end

γ = 1.0

Δt_values = 10.0 .^ LinRange(log(Δt_min) / log(10), log(Δt_max) / log(10), 16)

@threads for Δt in Δt_values
    constant_diffusion_compute_transition_times(scheme, Δt, n_transitions, set_tolerance; γ=γ)
end