module MatrixProductStates

# Basic functionality for contracting tensors in an efficient manner
include("tensorfunctions.jl")
export contract_tensors

# Functionality to provide some basic spin Hamiltonians as MPOs
include("operators.jl")
export getPauliMatrices, getIsingMPO, getHeisenbergMPO, getTotalSpinMPO, getIdentityMPO

# Methods implementing the MPS algorithm with open boundary conditions
include("MPS_OBC.jl")
export random_mps_obc, calculate_overlap, expectation_value, gaugeMPS, gaugeMPS!, apply_operator, find_groundstate, approximate_mps

end # module
