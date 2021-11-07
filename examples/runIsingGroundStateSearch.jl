############################################################################
# Find the ground state of the Ising Hamitlonian and get the magnetization #
############################################################################
using MatrixProductStates

# Number of Spins
N = 40
# Bond Dimension
D = 20
# Physical dimension
d = 2
# Desired Accuracy
acc = 1E-10
# Maximum number of sweeps
max_sweeps = 10

# Coupling between nearest neigbours
J = 1.0
# Coupling to external field
lambda = 1.0
# Generate Ising Hamiltonian
H = getIsingMPO(N, J, lambda)
# Get the MPO for the total spin
Sz = getTotalSpinMPO(N)

E0, mps, num_of_sweeps = find_groundstate(H, D, d, acc, max_sweeps)

println("Groundstate energy:            ", E0)
println("Groundstate energy density:    ", E0/N)
println("Number of sweeps:              ", num_of_sweeps)
println("Total spin:                    ", expectation_value(mps, Sz))