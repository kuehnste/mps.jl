##########################################################################
#       Some useful models and observables to get started with           #
##########################################################################

using LinearAlgebra

"""
    getPauliMatrices()

Provide the Pauli matrices as dense matrices.
"""
function getPauliMatrices()
    Id = [1.0 0;0.0 1.0]
    X = [0.0 1.0;1.0 0.0]
    Y = [0.0 -1.0im; 1.0im 0.0]
    Z = [1.0 0.0;0.0 -1.0]
    return Id, X, Y, Z
end


"""
    getHeisenbergMPO(N::Int64,J::Float64,lambda::Float64)
    
MPO implementation of the Heisenberg Hamiltonian 
H = J * sum_{i) (X^{i} X^{i+1} +Y^{i} Y^{i+1} + Z^{i} Z^{i+1}) + λ * sum_i X^i
where J, λ have to be real numbers .
"""
function getHeisenbergMPO(N::Int, J::Real, lambda::Real)::MPO{ComplexF64}
    # Provide the Pauli matrices
    Id, X, Y, Z = getPauliMatrices()
    # Initialize single tensor as complex array and the MPO as an array holding those
    first_tensor = zeros(ComplexF64, 1, 5, 2, 2)
    tensor = zeros(ComplexF64, 5, 5, 2, 2)
    last_tensor = zeros(ComplexF64, 5, 1, 2, 2)
    HeisenbergMPO = MPO{ComplexF64}(undef, N)

    # First MPO tensor
    first_tensor[1,1,:,:] = -lambda * Z
    first_tensor[1,2,:,:] = -J * X
    first_tensor[1,3,:,:] = -J * Y
    first_tensor[1,4,:,:] = -J * Z
    first_tensor[1,5,:,:] = Id

    # Tensors in the center of the MPO
    tensor[1,1,:,:] = Id
    tensor[2,1,:,:] = X
    tensor[3,1,:,:] = Y
    tensor[4,1,:,:] = Z
    tensor[5,1,:,:] = -lambda * Z
    tensor[5,2,:,:] = -J * X
    tensor[5,3,:,:] = -J * Y
    tensor[5,4,:,:] = -J * Z
    tensor[5,5,:,:] = Id

    # Last MPO tensor
    last_tensor[1,1,:,:] = Id
    last_tensor[2,1,:,:] = X
    last_tensor[3,1,:,:] = Y
    last_tensor[4,1,:,:] = Z
    last_tensor[5,1,:,:] = -lambda * Z

    # Fill the MPO
    HeisenbergMPO[1] = first_tensor
    for i = 2:N - 1
        HeisenbergMPO[i] = tensor
    end
    HeisenbergMPO[N] = last_tensor

    return HeisenbergMPO
end


"""
    getIsingMPO(N::Int64,J::Float64,lambda::Float64)
    
Function to generate the MPO for the Ising Hamilton of the form -J * sum_{i)X^i X^i+1 - lambda * sum_i Z^i.
"""
function getIsingMPO(N::Int, J::Real, lambda::Real)::MPO{Float64}
    # Provide the Pauli matrices
    Id, X, _, Z = getPauliMatrices()
    # Initialize single tensor as complex array and the MPO as an array holding those
    first_tensor = zeros(Float64, 1, 5, 2, 2)
    tensor = zeros(Float64, 5, 5, 2, 2)
    last_tensor = zeros(Float64, 5, 1, 2, 2)
    IsingMPO = MPO{Float64}(undef, N)

    # First MPO tensor
    first_tensor[1,1,:,:] = Id
    first_tensor[1,2,:,:] = -J * X
    first_tensor[1,3,:,:] = -lambda * Z

    # Tensors in the center of the MPO
    tensor[1,1,:,:] = Id
    tensor[1,2,:,:] = -J * X
    tensor[1,3,:,:] = -lambda * Z
    tensor[2,3,:,:] = X
    tensor[3,3,:,:] = Id

    # Last MPO tensor
    last_tensor[1,1,:,:] = -lambda * Z
    last_tensor[2,1,:,:] = X
    last_tensor[3,1,:,:] = Id

    # Fill the MPO
    IsingMPO[1] = first_tensor
    for i = 2:N - 1
        IsingMPO[i] = tensor
    end
    IsingMPO[N] = last_tensor

    return IsingMPO
end


"""
    getTotalSpinMPO(N::Int)
    
Function to generate the MPO for sum_{i)Z^i.
"""
function getTotalSpinMPO(N::Int)::MPO{Float64}
    # Provide the Pauli matrices
    Id, X, Y, Z = getPauliMatrices()
    # Initialize the result
    TotalSpinMPO = MPO{Float64}(undef, N)

    # Initialize single tensor as complex array and the MPO as an array holding those
    first_tensor = zeros(Float64, 1, 2, 2, 2)
    tensor = zeros(Float64, 2, 2, 2, 2)
    last_tensor = zeros(Float64, 2, 1, 2, 2)

    # First MPO tensor
    first_tensor[1,1,:,:] = Id
    first_tensor[1,2,:,:] = Z

    # Tensors in the center of the MPO
    tensor[1,1,:,:] = Id
    tensor[1,2,:,:] = Z
    tensor[2,2,:,:] = Id

    # Last MPO tensor
    last_tensor[1,1,:,:] = Z
    last_tensor[2,1,:,:] = Id

    # Fill the MPO
    TotalSpinMPO[1] = first_tensor
    for i = 2:N - 1
        TotalSpinMPO[i] = tensor
    end
    TotalSpinMPO[N] = last_tensor

    return TotalSpinMPO
end


"""
    getIdentityMPO(N::Int,d::Int,tensortype::Type{T}=Float64) where T
    
Function to generate the MPO for the Identity which is useful to check contractions and other functions.
"""
function getIdentityMPO(N::Int, d::Int, tensortype::Type{T}=Float64)::MPO{T} where T
    # Initialize single tensor as complex array and the MPO as an array holding those
    tensor = reshape(Matrix{T}(I, d, d), (1, 1, d, d))
    IdentityMPO = MPO{T}(undef, N)

    # Fill the MPO
    for i = 1:N
        IdentityMPO[i] = tensor
    end

    return IdentityMPO
end


