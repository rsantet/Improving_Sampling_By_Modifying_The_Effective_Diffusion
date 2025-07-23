using StaticArrays
using JLD2
using LinearAlgebra

# ========= System ========= #
n_particles_per_dimension::Int = 4
n_particles::Int = n_particles_per_dimension^2
d::Int = 2 * n_particles
σ::Float64 = 1.0 # R in the paper
r0::Float64 = 2^(1 // 6) * σ
density::Float64 = 0.7
box_length::Float64 = sqrt(n_particles / density)
cell_length::Float64 = box_length / n_particles_per_dimension
@assert r0 < box_length / 2 # check that a particle cannot interact with its periodized copy
β::Float64 = 1.0

# ========= Potentials ========= #

ε::Float64 = 1.0
w::Float64 = 0.35
z_tol::Float64 = 2.0 # Increase if CV values during test runs exceed this value
r1::Float64 = (box_length - 4 * w) / 4 # if one wants V_DW(0)=V_DW(box_length/2)...
@assert r1 + 2 * w * z_tol < box_length / 2 # check that the interactions due to the double-well potential let the dimer go to almost all values of the CV
# i.e. the potential has to be C^2 over the simulation...
h::Float64 = 2.0
σ6::Float64 = σ^6

function WCA(r)
    if r > r0
        return 0.0
    else
        r6 = r^6
        return 4 * ε * σ6 * (σ6 / r6 - 1) / r6 + ε
    end
end

function ∇WCA(r)
    if r > r0
        return 0.0
    else
        r6 = r^6
        return -24.0 * ε * σ6 * (2 * σ6 / r6 - 1) / r6 / r
    end
end

function DW(r)
    d = (r - r1 - w) / w
    return h * (1 - d^2)^2
end

function ∇DW(r)
    d = (r - r1 - w) / w
    return -4 * h * (1 - d^2) * d / w
end

function V(q, unit_vector, distance)
    tmp = 0.0

    # Solvent/Solvent interaction
    @inbounds for i in 3:n_particles
        @inbounds for j in i+1:n_particles
            length_PBC!(q, unit_vector, distance, j, i)
            tmp += WCA(distance[1])
        end # for j; solvent/solvent interaction
    end # for i; solvent/solvent interaction

    # Solvent/Dimer interaction
    @inbounds for i in 1:2
        @inbounds for j in 3:n_particles
            length_PBC!(q, unit_vector, distance, j, i)
            # println("SOLVENT/DIMER ", i, " ", j)
            tmp += WCA(distance[1])
        end # for j; solvent/Dimer interaction
    end # for i; solvent/Dimer interaction

    # Dimer/Dimer interaction
    length_PBC!(q, unit_vector, distance, 2, 1)
    tmp += DW(distance[1])

    return tmp
end

function reset_forces!(F)

    @inbounds for i in eachindex(F)
        F[i] = SVector(0.0, 0.0)
    end

end

function length_PBC!(q, unit_vector, distance, j, i)
    unit_vector[1] = q[j] - q[i]
    unit_vector[1] -= box_length * round.(unit_vector[1] / box_length)
    distance[1] = norm(unit_vector[1])
    unit_vector[1] /= distance[1]
end

function compute_forces!(q, F, unit_vector, distance)

    reset_forces!(F)

    # Solvent/Solvent interaction
    @inbounds for i in 3:n_particles
        @inbounds for j in i+1:n_particles
            length_PBC!(q, unit_vector, distance, j, i) # qj-qi
            act = ∇WCA(distance[1])
            F[i] += act * unit_vector[1]
            F[j] -= act * unit_vector[1]
        end # for j; solvent/solvent interaction
    end # for i; solvent/solvent interaction

    # Solvent/Dimer interaction
    @inbounds for i in 1:2
        @inbounds for j in 3:n_particles
            length_PBC!(q, unit_vector, distance, j, i)
            act = ∇WCA(distance[1])
            F[i] += act * unit_vector[1]
            F[j] -= act * unit_vector[1]
        end # for j; solvent/Dimer interaction
    end # for i; solvent/Dimer interaction

    # Dimer/Dimer interaction
    length_PBC!(q, unit_vector, distance, 2, 1) # q2-q1
    act = ∇DW(distance[1])
    F[1] += act * unit_vector[1]
    F[2] -= act * unit_vector[1]
end

function periodize!(q)
    for i in eachindex(q)
        q[i] -= box_length * floor.(q[i] / box_length)
    end
end

# ========= Collective Variable ========= #

function ξ(q, unit_vector, distance)
    length_PBC!(q, unit_vector, distance, 1, 2)
    return (distance[1] - r1) / 2 / w
end

function ∇ξ!(q, ∇z, unit_vector, distance)
    length_PBC!(q, unit_vector, distance, 1, 2)
    ∇z[1] = unit_vector[1] / 2 / w
    ∇z[2] = -∇z[1]
end

square_norm_∇z::Float64 = 1 / 2 / w^2

function ∇2ξ!(q, ∇2z, unit_vector, distance)
    length_PBC!(q, unit_vector, distance, 1, 2)

    ∇2z[1, 1] = 1 - unit_vector[1][1]^2
    ∇2z[1, 2] = -unit_vector[1][1] * unit_vector[1][2]
    ∇2z[1, 3] = -∇2z[1, 1]
    ∇2z[1, 4] = -∇2z[1, 2]

    ∇2z[2, 1] = ∇2z[1, 2]
    ∇2z[2, 2] = 1 - unit_vector[1][2]^2
    ∇2z[2, 3] = -∇2z[1, 2]
    ∇2z[2, 4] = -∇2z[2, 2]

    ∇2z[3, 1] = -∇2z[1, 1]
    ∇2z[3, 2] = -∇2z[1, 2]
    ∇2z[3, 3] = ∇2z[1, 1]
    ∇2z[3, 4] = ∇2z[1, 2]

    ∇2z[4, 1] = -∇2z[1, 2]
    ∇2z[4, 2] = -∇2z[2, 2]
    ∇2z[4, 3] = ∇2z[1, 2]
    ∇2z[4, 4] = ∇2z[2, 2]

    ∇2z ./= (2 * w * distance[1])

end

function Δξ!(q, unit_vector, distance)
    length_PBC!(q, unit_vector, distance, 1, 2)
    return 1 / w / distance[1]
end

# ========= Free energy and Mean force ========= #

# The z_values are moved by Δz/2 from z_min, z_max, only considering interior points
z_min::Float64 = -0.5
z_max::Float64 = 1.5
nz::Int = 100
Δz::Float64 = (z_max - z_min) / nz

z_values_fp = "data/TI/z_values.jld2"
FE_fp = "data/TI/free_energy.jld2"
MF_fp = "data/TI/mean_force.jld2"

if !(isfile(z_values_fp)& isfile(FE_fp) & isfile(MF_fp))
    throw("You should run thermodynamic integration (TI) first. If performing only adaptives schemes, comment the lines loading FE/MF objects and constructing the FM/MF vectors in input.jl.")
end

z_values::Vector{Float64} = load_object(z_values_fp)
FE::Vector{Float64} = load_object(FE_fp)
MF::Vector{Float64} = load_object(MF_fp)


# For stability reasons, we consider only a portion of the interval [z_min,z_max] for the learned Free Energy
z_min_FE = -0.25
z_max_FE = 1.25
index_z_min_FE = floor(Int, (z_min_FE - z_min) / Δz)
index_z_max_FE = floor(Int, (z_max_FE - z_min) / Δz)

FE = FE[index_z_min_FE:index_z_max_FE]
MF = MF[index_z_min_FE:index_z_max_FE]

# Ensuring that minimum(F) over [z_min_FE,z_max_FE] is 0
FE .-= minimum(FE)

compact_state = -0.05930930930930933
stretched_state = 0.9795795795795795

function mean_force(z)
    if z < z_min_FE || z > z_max_FE
        return 0.0
    end

    index = floor(Int, (z - z_min_FE) / Δz) + 1

    return MF[index]
end

function free_energy(z)
    if z <= z_min_FE
        return FE[1]
    elseif z >= z_max_FE
        return FE[end]
    else
        index = floor(Int, (z - z_min_FE) / Δz) + 1
        left_value = FE[index]

        left_z = z_min_FE + (index - 1) * Δz
        return left_value + MF[index] * (z - left_z)
    end
end


# ========= Diffusion ========= #

effective_diffusion_squared::Float64 = 1 / 2 / w^2

function κ(α)
    tmp = 0.0
    for i in eachindex(FE[1:end-1])
        tmp += sqrt(d - 1 + exp(2 * α * β * FE[i]) / effective_diffusion_squared^2) * exp(-β * FE[i]) * Δz
    end
    return 1 / tmp
end

function a(z, α)
    return exp(α * β * free_energy(z)) / effective_diffusion_squared
end

function ∇a(z, α)
    return α * β * mean_force(z) * a(z, α)
end

# Useful for MALA update
function compute_sqrt_div_inv_diffusion!(q, ∇z, projection_matrix, D, sqrt_D, div_D, inv_D, unit_vector, distance, α, κ_α)

    z = ξ(q, unit_vector, distance)
    ∇ξ!(q, ∇z, unit_vector, distance)

    a_z = a(z, α)
    ∇a_z = ∇a(z, α)

    projection_matrix .= outer_product(∇z, ∇z) / square_norm_∇z
    D .= κ_α * (I + (a_z - 1) * projection_matrix)
    sqrt_D .= sqrt(κ_α) * (I + (sqrt(a_z) - 1) * projection_matrix)
    inv_D .= (I + (1 / a_z - 1) * projection_matrix) / κ_α

    div_D .= κ_α * ((a_z - 1) * Δξ!(q, unit_vector, distance) / square_norm_∇z + ∇a_z) * ∇z

end

# ========= Adaptive Algorithm ========= #
n_minimum_number_of_samples::Int = 100
n_mean_force_update::Int = 1

function get_index_RC(z)
    if z < z_min || z > z_max
        return 0
    end
    return floor(Int, (z - z_min) / Δz) + 1
end

function get_mean_force(index::Int, MeanForce::Matrix)

    if index == 0 # Value of the RC do not lie in the considered range
        return 0.0
    end

    if MeanForce[index, 1] > n_minimum_number_of_samples # If the value has been seen sufficiently
        return MeanForce[index, 2] / MeanForce[index, 1]
    end

    return 0.0
end

function get_free_energy(index::Int, z, FreeEnergy::Vector, MeanForce::Matrix)

    if index == 0 # Value of the RC do not lie in the considered range
        return 0.0
    end

    # Linear interpolation
    left_value = FreeEnergy[index]
    left_z = z_min + (index - 1) * Δz
    slope = get_mean_force(index, MeanForce)

    return left_value + slope * (z - left_z)
end

function update_free_energy!(MeanForce::Matrix, FreeEnergy::Vector)

    # Reset the Free Energy values
    FreeEnergy .= 0

    # Simple quadrature formula
    for i in axes(MeanForce, 1)[1:end-1]
        if MeanForce[i, 1] > n_minimum_number_of_samples # If the value has been seen sufficiently
            FreeEnergy[i+1] = FreeEnergy[i] + MeanForce[i, 2] * Δz / MeanForce[i, 1]
        else
            FreeEnergy[i+1] = FreeEnergy[i]
        end
    end

    # Translation of the free energy so that is positive
    FreeEnergy .-= minimum(FreeEnergy)
end

function update_mean_force!(MeanForceUpdate::Matrix, MeanForce::Matrix)

    for i in axes(MeanForce, 1)
        MeanForce[i, 1] += MeanForceUpdate[i, 1]
        MeanForce[i, 2] += MeanForceUpdate[i, 2]
    end

    # Reset the update
    MeanForceUpdate .= 0
end

function update_κ_α(FreeEnergy::Vector, α)
    tmp = 0.0
    for i in eachindex(FreeEnergy[1:end-1])
        tmp += sqrt(d - 1 + exp(2 * α * β * FreeEnergy[i]) / effective_diffusion_squared^2) * exp(-β * FreeEnergy[i]) * Δz
    end
    κ_α = 1 / tmp
    sqrt_κ_α = sqrt(κ_α)
    return κ_α, sqrt_κ_α
end

function adaptive_a(FE, α)
    return exp(α * β * FE) / effective_diffusion_squared
end

function ∇adaptive_a(MF, FE, α)
    return α * β * MF * adaptive_a(FE, α)
end


function adaptive_compute_sqrt_div_inv_diffusion!(q, ∇z, projection_matrix, D, sqrt_D, div_D, inv_D, unit_vector, distance, MF, FE, α, κ_α, sqrt_κ_α)

    ∇ξ!(q, ∇z, unit_vector, distance)

    a_z = adaptive_a(FE, α)
    ∇a_z = ∇adaptive_a(MF, FE, α)

    projection_matrix .= outer_product(∇z, ∇z) / square_norm_∇z
    D .= κ_α * (I + (a_z - 1) * projection_matrix)
    sqrt_D .= sqrt_κ_α * (I + (sqrt(a_z) - 1) * projection_matrix)
    inv_D .= (I + (1 / a_z - 1) * projection_matrix) / κ_α

    div_D .= κ_α * ((a_z - 1) * Δξ!(q, unit_vector, distance) / square_norm_∇z + ∇a_z) * ∇z

end

# ========= Newton Algorithm ========= #
max_iterations_newton::Int = 1000
tol_cauchy_newton::Float64 = 1e-12
tol_root_newton::Float64 = 1e-12

# ========= Reversibility tolerance ========= #
tol_reversibility::Float64 = 1e-6

# ========= Transition time computations ========= #
set_tolerance = 0.1

function is_in_compact_form(
    q,
    unit_vector, distance, set_tolerance
)
    z = ξ(q, unit_vector, distance)
    return z < compact_state + set_tolerance
end

function is_in_expanded_form(
    q,
    unit_vector, distance, set_tolerance
)
    z = ξ(q, unit_vector, distance)
    return z > stretched_state - set_tolerance
end

# ========= Initializers ========= #
function initialize_OL(z, α) # for MALA
    q = reshape([SVector((0.5 + j) * cell_length, (0.5 + i) * cell_length) for i = 0:n_particles_per_dimension-1, j = 0:n_particles_per_dimension-1], n_particles)

    # Set dimer bond length such that CV value is z
    q[2] = SVector(q[1][1], q[1][2] + r1 + 2 * w * z)

    # Allocate structures for force computations
    F = zero(q)
    unit_vector = [SVector(0.0, 0.0)]
    distance = [0.0]
    compute_forces!(q, F, unit_vector, distance)

    # Allocate structures for collective variable
    ∇z = [SVector(0.0, 0.0) for i = 1:2]
    ∇ξ!(q, ∇z, unit_vector, distance)

    # Allocate structures for diffusion operator
    κ_α = κ(α)
    projection_matrix = Matrix(1.0I, 4, 4)
    D = Matrix(1.0I, 4, 4)
    sqrt_D = Matrix(1.0I, 4, 4)
    inv_D = Matrix(1.0I, 4, 4)
    div_D = [SVector(0.0, 0.0), SVector(0.0, 0.0)]

    compute_sqrt_div_inv_diffusion!(q, ∇z, projection_matrix, D, sqrt_D, div_D, inv_D, unit_vector, distance, α, κ_α)

    return q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, inv_D, div_D, κ_α
end

function initialize_adaptive_OL(z, α) # for adaptive MALA
    q = reshape([SVector((0.5 + j) * cell_length, (0.5 + i) * cell_length) for i = 0:n_particles_per_dimension-1, j = 0:n_particles_per_dimension-1], n_particles)

    # Set dimer bond length such that CV value is z
    q[2] = SVector(q[1][1], q[1][2] + r1 + 2 * w * z)

    # Allocate structures for force computations
    F = zero(q)
    unit_vector = [SVector(0.0, 0.0)]
    distance = [0.0]
    compute_forces!(q, F, unit_vector, distance)

    # Allocate structures for collective variable
    ∇z = [SVector(0.0, 0.0) for i = 1:2]
    ∇ξ!(q, ∇z, unit_vector, distance)

    # Allocate structures for diffusion operator
    projection_matrix = Matrix(1.0I, 4, 4)
    D = Matrix(1.0I, 4, 4)
    sqrt_D = Matrix(1.0I, 4, 4)
    inv_D = Matrix(1.0I, 4, 4)
    div_D = [SVector(0.0, 0.0), SVector(0.0, 0.0)]

    # Structure to hold the mean force over time
    # First axis is the histogram of the RCs
    # Second axis is the corresponding value of the mean force
    MeanForce = zeros(Float64, nz, 2)

    # Holds the values of the free energy
    FreeEnergy = zeros(Float64, nz)

    # Holds the update for the mean force
    MeanForceUpdate = zeros(Float64, nz, 2)

    κ_α, sqrt_κ_α = update_κ_α(FreeEnergy, α)

    z = ξ(q, unit_vector, distance)
    index = get_index_RC(z)
    # very first update
    lmf = (dot(F[2] - F[1], unit_vector[1]) - 2 / β / distance[1]) * w
    if index != 0
        MeanForceUpdate[index, 1] += 1
        MeanForceUpdate[index, 2] += lmf
    end

    MF = get_mean_force(index, MeanForce)
    FE = get_free_energy(index, z, FreeEnergy, MeanForce)

    adaptive_compute_sqrt_div_inv_diffusion!(q, ∇z, projection_matrix, D, sqrt_D, div_D, inv_D, unit_vector, distance, MF, FE, α, κ_α, sqrt_κ_α)

    return q, F, ∇z, unit_vector, distance, projection_matrix, D, sqrt_D, inv_D, div_D, MeanForce, FreeEnergy, MeanForceUpdate, κ_α, sqrt_κ_α
end

function initialize_HMC(z, α) # for RMHMC
    # Set initial configuration
    q = reshape([SVector((0.5 + j) * cell_length, (0.5 + i) * cell_length) for i = 0:n_particles_per_dimension-1, j = 0:n_particles_per_dimension-1], n_particles)

    # Set dimer bond length
    q[2] = SVector(q[1][1], q[1][2] + r1 + 2 * w * z)

    # Allocate structure for force computation
    F = zero(q)
    unit_vector = [SVector(0.0, 0.0)]
    distance = [0.0]

    compute_forces!(q, F, unit_vector, distance)

    # Allocate structures for collective variable
    ∇z = [SVector(0.0, 0.0) for i = 1:2]
    ∇2z = Matrix(1.0I, 4, 4)

    # Compute collective variable
    ∇ξ!(q, ∇z, unit_vector, distance)
    ∇2ξ!(q, ∇2z, unit_vector, distance)

    # Allocate structures for diffusion operator
    κ_α = κ(α)
    projection_matrix = Matrix(1.0I, 4, 4)
    inv_sqrt_D = Matrix(1.0I, 4, 4)

    # Compute diffusion matrices
    projection_matrix .= outer_product(∇z, ∇z) / square_norm_∇z
    inv_sqrt_D .= (I + (1 / sqrt(a(z, α)) - 1) * projection_matrix) / sqrt(κ_α)

    # Sample momenta
    p = randn(SVector{2,Float64}, n_particles) / sqrt(β)
    p[1:2] .= inv_sqrt_D * p[1:2]
    p[3:end] /= sqrt(κ_α)

    return q, p, F, ∇z, ∇2z, unit_vector, distance, projection_matrix, inv_sqrt_D, κ_α
end

# ========= Counter struct ========= #
mutable struct Counter
    MH_rejection::Int64
    forward_no_convergence_momenta_update::Int64
    forward_no_convergence_position_update::Int64
    backward_no_convergence_momenta_update::Int64
    backward_no_convergence_position_update::Int64
    no_reversibility::Int64
end

# ========= StaticArrays operations ========= #

import Base.*
function *(mat::Matrix, vec::Vector{SVector{2,T}}) where {T}
    flattened_vec = reduce(vcat, vec)
    @assert size(mat)[end] == length(flattened_vec)
    res = mat * flattened_vec
    return [SVector(res[i], res[i+1]) for i in 1:2:length(flattened_vec)]
end

import Base.\
function \(mat::Matrix, vec::Vector{SVector{2,T}}) where {T}
    flattened_vec = reduce(vcat, vec)
    @assert size(mat)[end] == length(flattened_vec)
    res = mat \ flattened_vec
    return [SVector(res[i], res[i+1]) for i in 1:2:length(flattened_vec)]
end

# Returns a matrix
function outer_product(X::Vector{SVector{2,T}}, Y::Vector{SVector{2,T}}) where {T}
    return reduce(vcat, X) * reduce(vcat, Y)'
end

import LinearAlgebra.dot
function dot(X::Vector{SVector{2,T}}, Y::Vector{SVector{2,T}}) where {T}
    return sum(map((x, y) -> x[1] * y[1] + x[2] * y[2], X, Y))
end

function square_norm(X::Vector{SVector{2,T}}) where {T}
    return dot(X, X)
end