############################################################################
# Find the ground state of the Ising Hamitlonian and get the magnetization #
############################################################################
using MatrixProductStates

let
    # Number of Spins
    N = 20
    # Bond Dimension
    D = 20
    # Physical dimension
    d = 2
    # Desired Accuracy
    acc = 1E-10
    # The number of time steps we take
    nsteps = 200
    # The size of each time step
    dt = 5E-2

    # Coupling between nearest neigbours
    J = 1.0
    # Coupling to external field
    lambda = 1.0
    # Generate Ising Hamiltonian H = -J * sum_{i=1}^{N-1} X^i X^i+1 - lambda * sum_{i=1}^N Z^i.
    H = getIsingMPO(N, J, lambda)
    # Get the MPO for the total spin
    Sz = getTotalSpinMPO(N)

    #=
    In the following we prepare the time evolution operators in MPO form using an
    odd-even decomposition for the evolution operator. Since the Ising Hamiltonian is
    of the form H = \sum_i h_{i,i+1} we can approximate the matrix exponential using 
    a Suzuki-Trotter decomposition as follows
    exp(-dt * H) ≈ \prod_{i even} exp(-dt * h_{i,i+1}) \prod_{i odd} exp(-dt * h_{i,i+1}).
    Since all terms starting at an odd (even) site commute, we can find an MPO
    representation by decomping one of the local terms exp(-dt * h_{i,i+1}) into an MPO,
    and using the same tensors at every site. Special care has to be taken on the
    boundaries. 
    =#

    # The required Pauli matrices for building the evolution MPOs
    Id, X, _, Z = getPauliMatrices()
    # Prepare the matrix representation for the different exponentials appearing in the odd-even decomposition
    Uloc_left_boundary = exp(-dt * (-J * kron(X, X) - lambda * kron(Z, Id) - 0.5 * lambda * kron(Id, Z)))
    Uloc_center = exp(-dt * (-J * kron(X, X) - 0.5 * lambda * (kron(Z, Id) + kron(Id, Z))))
    Uloc_right_boundary = exp(-dt * (-J * kron(X, X) - 0.5 * lambda * kron(Z, Id) - lambda * kron(Id, Z)))
    # The identity operator as MPO tensor
    id_tensor = reshape(Id, (1, 1, 2, 2))
    # Find the mpo representation for the local terms
    mpo_Uloc_left_boundary = decompose_into_mpo(Uloc_left_boundary, d)
    mpo_Uloc_center = decompose_into_mpo(Uloc_center, d)
    mpo_Uloc_right_boundary = decompose_into_mpo(Uloc_right_boundary, d)
    # Build the MPO for evolving the odd sites
    Uodd = MPO{Float64}(undef, N)
    for i = 1:2:N-1
        if i == 1
            Uodd[i] = mpo_Uloc_left_boundary[1]
            Uodd[i+1] = mpo_Uloc_left_boundary[2]
        elseif i == N - 1
            Uodd[i] = mpo_Uloc_right_boundary[1]
            Uodd[i+1] = mpo_Uloc_right_boundary[2]
        else
            Uodd[i] = mpo_Uloc_center[1]
            Uodd[i+1] = mpo_Uloc_center[2]
        end
    end
    if isodd(N)
        Uodd[N] = id_tensor
    end
    # Build the MPO for evolving the even sites
    Ueven = MPO{Float64}(undef, N)
    Ueven[1] = id_tensor
    for i = 2:2:N-1
        if i == N - 1
            Ueven[i] = mpo_Uloc_right_boundary[1]
            Ueven[i+1] = mpo_Uloc_right_boundary[2]
        else
            Ueven[i] = mpo_Uloc_center[1]
            Ueven[i+1] = mpo_Uloc_center[2]
        end
    end
    if iseven(N)
        Ueven[N] = id_tensor
    end

    # Evolve in imaginary time starting from a random, normalized state
    psi = random_mps_obc(N, D, d, Float64)
    gaugeMPS!(psi, :left, true)
    for i = 1:nsteps
        print("\rStep: ", i)
        # Apply the evolution operators
        psi = apply_operator(Uodd, psi)
        psi = apply_operator(Ueven, psi)
        # Compress and renormalize
        psi = svd_compress_mps(psi, D)
        gaugeMPS!(psi, :left, true)
    end
    E0 = expectation_value(psi, H)

    # Print the results
    println(" ")
    println("Groundstate energy:            ", E0)
    println("Groundstate energy density:    ", E0 / N)
    println("Total spin:                    ", expectation_value(psi, Sz))
    println("Entropy along center bond      ", compute_entropy(psi))

    nothing
end