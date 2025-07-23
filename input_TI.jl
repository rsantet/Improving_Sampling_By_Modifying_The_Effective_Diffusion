using StaticArrays
using JLD2
using LinearAlgebra

# ========= Thermodynamic Integration ========= #
TI_tolerance::Float64 = 1e-6

function local_mean_force(q, F, unit_vector, distance)
    length_PBC!(q, unit_vector, distance, 1, 2)
    return (dot(F[2] - F[1], unit_vector[1]) - 2 / β / distance[1]) * w
end

# ========= System ========= #
n_particles_per_dimension::Int = 4
n_particles::Int = n_particles_per_dimension^2
d::Int = 2 * n_particles
σ::Float64 = 1.0
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
# i.e. the potential has to be C^1 (and even more) over the simulation...
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

function initialize_OL(z) # for MALA
    q = reshape([SVector((0.5 + j) * cell_length, (0.5 + i) * cell_length) for i = 0:n_particles_per_dimension-1, j = 0:n_particles_per_dimension-1], n_particles)

    # Set dimer bond length such that CV value is z
    q[2] = SVector(q[1][1], q[1][2] + r1 + 2 * w * z)

    # Allocate structures for force computations
    F = zero(q)
    unit_vector = [SVector(0.0, 0.0)]
    distance = [0.0]

    # Compute forces
    compute_forces!(q, F, unit_vector, distance)

    return q, F, unit_vector, distance
end