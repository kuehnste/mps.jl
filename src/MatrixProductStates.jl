module MatrixProductStates

# Basic functionality for contracting tensors in an efficient manner
include("tensorfunctions.jl")
export contract_tensors

# Functionality to provide some basic spin Hamiltonians as MPOs
include("operators.jl")
export getPauliMatrices, getIsingMPO, getHeisenbergMPO, getTotalSpinMPO, getIdentityMPO

# Methods implementing the MPS algorithm with open boundary conditions
include("MPS_OBC.jl")
export Site, MPS, Operator, MPO, random_mps_obc, basis_state_obc, calculate_overlap, expectation_value, gaugeMPS, gaugeMPS!, apply_operator, apply_operator!, sum_states, sum_operators, contract_virtual_indices, find_groundstate, approximate_mps, compute_entropy, sample_from_mps!, sample_from_mps, svd_compress_mps

end # module
