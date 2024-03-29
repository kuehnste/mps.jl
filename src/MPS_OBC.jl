##########################################################################
#   	Basic functionality for doing MPS calculations with OBC     	#
##########################################################################
using LinearAlgebra
using Arpack

# Define sites as 3 legged tensors filled with elements of some type
Site{T} = Array{T,3}
# An MPS is an array of Sites of a certain type
MPS{T} = Vector{Site{T}}
# Define operators as 4 legged tensors filled with elements of some type
Operator{T} = Array{T,4}
# A MPO is an array of Operators
MPO{T} = Vector{Operator{T}}

"""
    random_mps_obc(N::Int, D::Int, d, tensortype::Type{T}=ComplexF64)::MPS{T} where T
    
Function to generate a random MPS with OBC
N: Number of sites
D: Bond dimension
d: Vector of physical dimensions, if a single integer dimension is given, it is assumed that all physical dimensions are the same, if varying dimensions are wanted a vector of integers of length N has to be supplied
The convention for the 3 index tensors is that the first two are the virtual ones, the last one is the physical one
   3
   |
1--A--2
"""
function random_mps_obc(N::Int, D::Int, d, tensortype::Type{T} = ComplexF64)::MPS{T} where {T}
    # Ensure the the system size is at least two
    @assert(N > 1)

    # Initialize
    mps = Array{Site{T}}(undef, N)

    if isa(d, Number)
        # Single number, I assume all dimensions are the same
        dim = d * ones(Int64, N)
    else
        # Make sure that the input is really a vector as julia distinguishes between nx1 matrices and column vectors
        dim = vec(d)
        # Make sure input is valid
        @assert(length(dim) == N)
    end

    # Left boundary tensor (row vector)
    mps[1] = rand(tensortype, 1, D, dim[1])

    # Right boundary tensor (column vector)
    mps[N] = rand(tensortype, D, 1, dim[N])

    # Tensors in between
    for i = 2:N-1
        mps[i] = rand(tensortype, D, D, dim[i])
    end

    return mps
end

"""
    basis_state_obc(configuration::Vector{<:Int}, d::Int=2)::MPS

Prepare the a product state corresponding |configuration> on N sites where configuration is an Array containing N elements from 1 to d. Each site is then initialized in the dth canonical basis state
"""
function basis_state_obc(configuration::Vector{<:Int}, d::Int = 2)::MPS{Float64}
    # Some error checking
    if any(x -> (x < 1 || x > d), configuration)
        throw(ArgumentError("configuration must contain integer elements in the range from 1 to d, got d=$(repr(d)), configuration=$(repr(configuration))"))
    end
    # Generate the MPS
    N = length(configuration)
    psi = MPS{Float64}(undef, N)
    tensors = Vector{Array{Float64,3}}(undef, d)
    for i = 1:d
        tmp = zeros(Float64, 1, 1, d)
        tmp[1, 1, i] = 1.0
        tensors[i] = tmp
    end
    for i = 1:N
        psi[i] = tensors[configuration[i]]
    end
    return psi
end

"""
    calculate_overlap(mps1::MPS,mps2::MPS)::Number

Given two mps, compute the overlap <mps1|mps2>.
"""
function calculate_overlap(mps1::MPS, mps2::MPS)::Number
    # Make sure, that the two MPS have the same length
    N1 = length(mps1)
    N2 = length(mps2)
    @assert(N1 == N2)

    # Now compute ther overlap
    overlap = ones(Float64, 1, 1)
    for i = 1:N1
        overlap = contract_tensors(overlap, [2], mps2[i], [1])
        overlap = contract_tensors(conj(mps1[i]), [1; 3], overlap, [1; 3])
    end

    return overlap[1]
end


"""
   expectation_value(mps::MPS, mpo::MPO)::Number
   
Given a MPS and a MPO compute the expectation value <mps|mpo|mps>.
"""
function expectation_value(mps::MPS, mpo::MPO)::Number
    # Make sure, that the two MPS have the same length
    N1 = length(mps)
    N2 = length(mpo)
    @assert(N1 == N2)

    # Contract zipper like
    val = ones(Float64, 1, 1, 1)
    for i = 1:N1
        val = contract_tensors(val, [1], conj(mps[i]), [1])
        val = contract_tensors(val, [1; 4], mpo[i], [1; 3])
        val = contract_tensors(val, [1; 4], mps[i], [1; 3])
    end

    return val[1]
end


"""
    gaugeMPS(mps::MPS{T}, direction::Symbol=:right, normalize::Bool=false)::MPS{T}
    
Function to bring an MPS with OBC in canonical form
direction: tells, whether it will be left or right canonical gauge
If normalize is set to true, the resulting state will be normalized.
"""
function gaugeMPS(mps::MPS{T}, direction::Symbol = :right, normalize::Bool = false)::MPS{T} where {T}
    mps_gauged = deepcopy(mps)
    gaugeMPS!(mps_gauged, direction, normalize)
    return mps_gauged
end


"""
    gaugeMPS!(mps::MPS,direction::Dir=right,)
    
Bring an MPS with OBC in canonical form. Direction tells, whether it will be left or right canonical gauge If normalize is set to true, the resulting state will be
normalized. This function overwrites the input MPS with its gauged version
"""
function gaugeMPS!(mps::MPS, direction::Symbol = :right, normalize::Bool = false)
    # Check that we got a meaningful direction
    if (direction != :left && direction != :right)
        throw(ArgumentError("direction must be :left or :right, got $(repr(direction))"))
    end
    # Start gauging
    N = length(mps)
    if (direction == :left)
        # Case that I want left canonical gauge
        M = mps[1]
        for i = 1:N-1
            mps[i], res = gauge_site(M, direction)
            M = contract_tensors(res, [2], mps[i+1], [1])
            if (i == (N - 1))
                if (normalize)
                    mps[N], _ = gauge_site(M, direction)
                else
                    mps[N] = M
                end
            end
        end
    else
        # Case that I want right canonical gauge
        M = mps[N]
        for i = N:-1:2
            mps[i], res = gauge_site(M, direction)
            M = contract_tensors(mps[i-1], [2], res, [1], [1; 3; 2])
            if (i == 2)
                if (normalize)
                    mps[1], _ = gauge_site(M, direction)
                else
                    mps[1] = M
                end
            end
        end
    end
end


"""
    gauge_site(A::Site{T}, direction::Dir)::Tuple{Site{T},Matrix{T}} where T

Bring a single tensor A of an MPS with OBC in canonical form
direction tells, whether it will be left or right normalized
The new tensor is M, the residual matrix which has to be multiplied 
in the next tensor is stored in res
"""
function gauge_site(A::Site{T}, direction::Symbol)::Tuple{Site{T},Matrix{T}} where {T}
    if (direction == :left)
        Dl, Dr, d = size(A)
        M = permutedims(A, [3 1 2])
        M = reshape(M, (d * Dl, Dr))
        U, S, V = svd(M)
        dsv = length(S)
        U = reshape(U, (d, Dl, dsv))
        M = permutedims(U, [2 3 1])
        res = Matrix(Diagonal(S)) * V'
        return M, res
    else
        Dl, Dr, d = size(A)
        M = permutedims(A, [1 3 2])
        M = reshape(M, (Dl, Dr * d))
        U, S, V = svd(M)
        dsv = length(S)
        V = V'
        V = reshape(V, (dsv, d, Dr))
        M = permutedims(V, [1 3 2])
        res = U * Matrix(Diagonal(S))
        return M, res
    end
end


"""
    contract_virtual_indices(mps::MPS)::Vector{<:Number}

Given an MPS contract the virtual indices such that one obtains a dense vector. The indices are ordered such that they are compatible with the standard Julia kronecker product.

Warning: the object constructed will have exponential memory requirements in terms of the number of sites, use with care!
"""
function contract_virtual_indices(mps::MPS)::Vector{<:Number}
    N = length(mps)

    # Since we deal with open boundary conditions, we drop the dummy indices one on the left (right) boundary for the first (last) tensor manually. We start from the right to have the physical indices in the order compatible with Julia's kronecker product
    res = mps[N][:, 1, :]
    for i = N-1:-1:2
        res = contract_tensors(res, [ndims(res) - 1], mps[i], [2])
    end
    res = contract_tensors(res, [ndims(res) - 1], mps[1][1, :, :], [1])

    # Now reshape the result accordingly
    res = reshape(res, prod(size(res)))

    return res
end


"""
    contract_virtual_indices(mps::MPO)::Matrix{<:Number}

Given an MPO contract the virtual indices such that one obtains a dense matrix. The indices are ordered such that they are compatible with the standard Julia kronecker product.

Warning: the object constructed will have exponential memory requirements in terms of the number of sites, use with care!
"""
function contract_virtual_indices(mpo::MPO)::Matrix{<:Number}
    N = length(mpo)

    # Since we deal with open boundary conditions, we drop the dummy indices one on the left (right) boundary for the first (last) tensor manually. We start from the right to have the physical indices in the order compatible with Julia's kronecker product
    res = mpo[N][:, 1, :, :]
    for i = N-1:-1:2
        res = contract_tensors(res, [ndims(res) - 2], mpo[i], [2])
    end
    res = contract_tensors(res, [ndims(res) - 2], mpo[1][1, :, :, :], [1])

    # Reshuffle the indices to be compatible with Julia's standard kronecker product
    res = permutedims(res, [collect(1:2:ndims(res)); collect(2:2:ndims(res))])

    # Now reshape the result accordingly
    dims = size(res)
    dr = prod(dims[1:N])
    dc = prod(dims[N+1:end])
    res = reshape(res, (dr, dc))

    return res
end


"""
    find_groundstate(H::Array{Operator},D::Int64,acc::Float64,max_sweeps::Int64=-1)

Given a Hamiltonian MPO H, find an MPS approximation for its ground state with bond 
dimension D converged to a relative accuracy acc. If a positive number for max_sweeps
is given, the maximum amount of iterations is limited to max_sweeps
"""
function find_groundstate(H::Vector{Operator{T}}, D::Int64, d::Int64, acc::Float64, max_sweeps::Int64) where {T}
    N = length(H)

    # Get the physical dimensions from the Hamiltonian
    d = Vector{Int64}(undef, length(H))
    for i = 1:length(H)
        d[i] = size(H[i], 3)
    end

    # Random MPS to start the calculation with
    mps = random_mps_obc(N, D, d)

    # Put it in right canonical gauge
    gaugeMPS!(mps, :right)
    # Precalculate the partial contractions
    LR = setup_R(H, mps)

    # Now start the sweeping
    num_of_sweeps = 0
    E0 = 0
    Eold = 1E5
    while (true)
        num_of_sweeps = num_of_sweeps + 1
        println("Sweep number ", num_of_sweeps)

        # From left to right starting by 1 up to N-1
        E_local = 0
        res = 0
        for i = 1:N-1
            Left = LR[i]
            Right = LR[i+1]
            E_local, M = solve_eigenvalue_problem(Left, Right, H[i])
            mps[i], _ = gauge_site(M, :left)
            LR[i+1] = update_left(Left, mps[i], mps[i], H[i])
        end

        # From left to right starting by 1 up to N-1
        for i = N:-1:2
            Left = LR[i]
            Right = LR[i+1]
            E_local, M = solve_eigenvalue_problem(Left, Right, H[i])
            mps[i], res = gauge_site(M, :right)
            LR[i] = update_right(Right, mps[i], mps[i], H[i])
        end

        # Check convergence criterion
        if (abs((Eold - E_local) / Eold) < acc)
            E0 = E_local
            # Restore normalisation (permutation is needed to bring the indices back into right order after contracting)
            mps[1] = contract_tensors(mps[1], [2], res, [1], [1; 3; 2])
            println("Convergence to desired accuracy achieved")
            break
        elseif (max_sweeps > 0 && num_of_sweeps >= max_sweeps)
            E0 = E_local
            # Restore normalisation (permutation is needed to bring the indices back into right order after contracting)
            mps[1] = contract_tensors(mps[1], [2], res, [1], [1; 3; 2])
            warn("Reached maximum number of iterations before convergence to desired accuracy")
            break
        end
        Eold = E_local
    end

    return E0, mps, num_of_sweeps
end


"""
    getHeff(Left,Right,W::Operator)
    
Construct the effective Hamiltonian
"""
function getHeff(Left, Right, W::Operator)
    Htemp = contract_tensors(W, [2], Right, [2])
    Htemp = contract_tensors(Left, [2], Htemp, [1])
    Htemp = permutedims(Htemp, [1; 5; 3; 2; 6; 4])
    dim1, dim2, dim3, dim4, dim5, dim6 = size(Htemp)
    Heff = reshape(Htemp, (dim1 * dim2 * dim3, dim4 * dim5 * dim6))

    return Heff
end


"""
    solve_eigenvalue_problem(Left,Right,W)
    
Function to construct and solve the eigenvalue problem which arises on each site for the effective Hamiltonian
"""
function solve_eigenvalue_problem(Left, Right, W)
    Htemp = contract_tensors(W, [2], Right, [2])
    Htemp = contract_tensors(Left, [2], Htemp, [1])
    Htemp = permutedims(Htemp, [1; 5; 3; 2; 6; 4])
    dim1, dim2, dim3, dim4, dim5, dim6 = size(Htemp)
    Heff = reshape(Htemp, (dim1 * dim2 * dim3, dim4 * dim5 * dim6))

    E_local, M = eigs(Heff, nev = 1, which = :SR)
    M = reshape(M, (dim1, dim5, dim3))

    return real(E_local[1]), M
end


"""
    setup_R(H::MPO{T1},mps::MPS{T2}) where {T1,T2}
    
Function to calculate the partial contractions needed to form the effective Hamiltonian
"""
function setup_R(H::MPO{T1}, mps::MPS{T2}) where {T1,T2}
    N = length(H)
    Tres = Base.return_types(*, (T1, T2))[1]
    LR = Vector{Array{Tres,3}}(undef, N + 1)

    # We need only N-1 partial contractions, however, we set the edges to dummy values 1 that we can recursively compute every contraction resuing the previous ones
    LR[1] = ones(Tres, 1, 1, 1)
    LR[N+1] = ones(Tres, 1, 1, 1)

    # Now compute the partial contractions starting from the right (as I start sweeping on the left in the ground state search)
    for i = N:-1:2
        tmp = contract_tensors(LR[i+1], [3], mps[i], [2])
        tmp = contract_tensors(H[i], [2; 4], tmp, [2; 4])
        LR[i] = contract_tensors(conj(mps[i]), [2; 3], tmp, [3; 2])
    end

    return LR
end


"""
    update_left(LR, M_left::Site, M_right::Site, W::Operator)
    
Function to perform an update the partial contractions required for iterative ground state search starting from the left stored in LR 
"""
function update_left(LR, M_left::Site, M_right::Site, W::Operator)
    res = contract_tensors(LR, [3], M_right, [1])
    res = contract_tensors(W, [1; 4], res, [2; 4])
    res = contract_tensors(conj(M_left), [1; 3], res, [3; 2])

    return res
end


"""
    update_right(LR, M_left::Site, M_right::Site, W::Operator)

Function to perform an update the partial contractions required for iterative ground state search starting from the right stored in LR 
"""
function update_right(LR, M_left::Site, M_right::Site, W::Operator)
    res = contract_tensors(LR, [3], M_right, [2])
    res = contract_tensors(W, [2; 4], res, [2; 4])
    res = contract_tensors(conj(M_left), [2; 3], res, [3; 2])

    return res
end


"""
    apply_operator(operator::MPO{T1}, mps::MPS{T2})::MPS where {T1,T2}
    
Apply an operator given as MPO to an MPS. The resulting MPS will have a bond dimension that is the product of the bond dimensions of the MPS and the MPO.
"""
function apply_operator(operator::MPO{T1}, mps::MPS{T2})::MPS where {T1,T2}
    N1 = length(mps)
    N2 = length(operator)
    @assert(N1 == N2)

    # Generate a new MPO of the correct type
    Tres = Base.return_types(*, (T1, T2))[1]
    res = MPS{Tres}(undef, N1)

    # Apply the MPO to the MPS and generate new MPS
    for i = 1:N1
        temp = contract_tensors(operator[i], [4], mps[i], [3])
        temp = permutedims(temp, (1, 4, 2, 5, 3))
        dim1, dim2, dim3, dim4, dim5 = size(temp)
        res[i] = reshape(temp, (dim1 * dim2, dim3 * dim4, dim5))
    end
    return res
end


"""
    apply_operator!(operator::MPO{T1}, mps::MPS{T2})::MPS where {T1,T2}
    
Apply an operator given as MPO to an MPS. The resulting MPS will have a bond dimension that is the product of the bond dimensions of the MPS and the MPO. The input will be overwritten by the result. This requires that the type of the elements of the input MPS is able to accomodate the result of multiplying the MPO tensors into the MPS tensors (e.g. applying a complex MPO to a real MPS cannot be done inpalce as the result will be complex)
"""
function apply_operator!(operator::MPO, mps::MPS)
    N1 = length(mps)
    N2 = length(operator)
    @assert(N1 == N2)

    # Contract the tensors for each site
    for i = 1:N1
        temp = contract_tensors(operator[i], [4], mps[i], [3])
        temp = permutedims(temp, (1, 4, 2, 5, 3))
        dim1, dim2, dim3, dim4, dim5 = size(temp)
        mps[i] = reshape(temp, (dim1 * dim2, dim3 * dim4, dim5))
    end
end


"""
    apply_operator(op1::MPO{T1}, op2::MPO{T2})::MPO where {T1,T2}
    
Multiply two MPOs together to get an expression for op2 * op1 in MPO form. The resulting MPS will have a bond dimension that is the product of the bond dimensions of both MPOs.
"""
function apply_operator(op1::MPO{T1}, op2::MPO{T2})::MPO where {T1,T2}
    N1 = length(op1)
    N2 = length(op2)
    @assert(N1 == N2)

    # Generate a new MPS of the correct type
    Tres = Base.return_types(*, (T1, T2))[1]
    res = MPO{Tres}(undef, N1)

    # Contract the tensors for each site
    for i = 1:N1
        temp = contract_tensors(op2[i], [4], op1[i], [3])
        temp = permutedims(temp, (1, 4, 2, 5, 3, 6))
        dim1, dim2, dim3, dim4, dim5, dim6 = size(temp)
        res[i] = reshape(temp, (dim1 * dim2, dim3 * dim4, dim5, dim6))
    end
    return res
end


"""
    sum_states(mps1::MPS{T1}, mps2::MPS{T2})::MPS where {T1,T2}
    
Add two MPSs together to get an expression for mps1 + mps2 in MPS form. The resulting MPS will have a bond dimension that is the sum of the bond dimensions of both MPOs.
"""
function sum_states(mps1::MPS{T1}, mps2::MPS{T2})::MPS where {T1,T2}
    N1 = length(mps1)
    N2 = length(mps2)
    @assert(N1 == N2)

    # Generate a new MPO of the correct type
    Tres = Base.return_types(+, (T1, T2))[1]
    res = MPS{Tres}(undef, N1)

    # The first site needs special treatment
    tensor1 = mps1[1]
    tensor2 = mps2[1]
    _, Dr1, d1 = size(tensor1)
    _, Dr2, d2 = size(tensor2)
    new_tensor = zeros(Tres, 1, Dr1 + Dr2, d1)
    for r = 1:d1
        new_tensor[1, :, r] = [tensor1[1:1, :, r] tensor2[1:1, :, r]]
    end
    res[1] = new_tensor
    # The tensors in between
    for i = 2:N1-1
        tensor1 = mps1[i]
        tensor2 = mps2[i]
        Dl1, Dr1, d1 = size(tensor1)
        Dl2, Dr2, d2 = size(tensor2)
        new_tensor = zeros(Tres, Dl1 + Dl2, Dr1 + Dr2, d1)
        for r = 1:d1
            new_tensor[1:Dl1, 1:Dr1, r] = tensor1[:, :, r]
            new_tensor[Dl1+1:end, Dr1+1:end, r] = tensor2[:, :, r]
        end
        res[i] = new_tensor
    end
    # The last site needs special treatment
    tensor1 = mps1[N1]
    tensor2 = mps2[N1]
    Dl1, _, d1 = size(tensor1)
    Dl2, _, d2 = size(tensor2)
    new_tensor = zeros(Tres, Dl1 + Dl2, 1, d1)
    for r = 1:d1
        new_tensor[:, 1, r] = [tensor1[:, 1, r]; tensor2[:, 1, r]]
    end
    res[N1] = new_tensor

    return res
end


"""
    sum_operators(op1::MPO{T1}, op2::MPO{T2})::MPO where {T1,T2}
    
Add two MPOs together to get an expression for op2 + op1 in MPO form. The resulting MPO will have a bond dimension that is the sum of the bond dimensions of both MPOs.
"""
function sum_operators(op1::MPO{T1}, op2::MPO{T2})::MPO where {T1,T2}
    N1 = length(op1)
    N2 = length(op2)
    @assert(N1 == N2)

    # Generate a new MPO of the correct type
    Tres = Base.return_types(+, (T1, T2))[1]
    res = MPO{Tres}(undef, N1)

    # The first site needs special treatment
    tensor1 = op1[1]
    tensor2 = op2[1]
    _, Dr1, dr1, dc1 = size(tensor1)
    _, Dr2, dr2, dc2 = size(tensor2)
    new_tensor = zeros(Tres, 1, Dr1 + Dr2, dr1, dc1)
    for r = 1:dr1
        for c = 1:dc1
            new_tensor[1, :, r, c] = [tensor1[1:1, :, r, c] tensor2[1:1, :, r, c]]
        end
    end
    res[1] = new_tensor
    # The tensors in between
    for i = 2:N1-1
        tensor1 = op1[i]
        tensor2 = op2[i]
        Dl1, Dr1, dr1, dc1 = size(tensor1)
        Dl2, Dr2, dr2, dc2 = size(tensor2)
        new_tensor = zeros(Tres, Dl1 + Dl2, Dr1 + Dr2, dr1, dc1)
        for r = 1:dr1
            for c = 1:dc1
                new_tensor[1:Dl1, 1:Dr1, r, c] = tensor1[:, :, r, c]
                new_tensor[Dl1+1:end, Dr1+1:end, r, c] = tensor2[:, :, r, c]
            end
        end
        res[i] = new_tensor
    end
    # The last site needs special treatment
    tensor1 = op1[N1]
    tensor2 = op2[N1]
    Dl1, _, dr1, dc1 = size(tensor1)
    Dl2, _, dr2, dc2 = size(tensor2)
    new_tensor = zeros(Tres, Dl1 + Dl2, 1, dr1, dc1)
    for r = 1:dr1
        for c = 1:dc1
            new_tensor[:, 1, r, c] = [tensor1[:, 1, r, c]; tensor2[:, 1, r, c]]
        end
    end
    res[N1] = new_tensor

    return res
end


"""
    compute_entropy(mps::MPS, n:Int)::Float64

Compute the von Neumann along for the bipartion of the sites into two subsets A={1,...,n} and B={n+1,...,N} where N is the length of the MPS. If n<=0 is supplied, a bipartion of into A. For the result to make sense, mps has to be a normalized quantum state.
"""
function compute_entropy(mps::MPS{T}, n::Int = 0)::Float64 where {T}
    # Extact the length and check if the given position is reasonable
    N = length(mps)
    @assert(n <= N)

    # The entropy of the entire state is simply zero since it is a pure state, so nothing to compute
    if n == N
        return 0.0
    end

    # In case a value n<=0 is given, we assume the bipartion is taken in the center
    if n <= 0
        n = Int(round(N / 2))
    end

    # Get a copy of the input
    mps_loc = deepcopy(mps)

    # Put the sites left to n into left canonical gauge    
    M = mps_loc[1]
    for i = 1:n
        mps_loc[i], res = gauge_site(M, :left)
        M = contract_tensors(res, [2], mps_loc[i+1], [1])
    end
    mps_loc[n+1] = M

    # Now start contracting from the right boundary to obtain the reduced density operator 
    rdm = ones(T, 1, 1)
    for i = N:-1:n+1
        rdm = contract_tensors(rdm, [2], mps_loc[i], [2])
        rdm = contract_tensors(conj(mps_loc[i]), [2; 3], rdm, [1; 3])
    end

    # Now diagonalize the reduced density matrix and compute the entropy    
    ev = real(eigvals(rdm))
    # Eigenvalues that are numerically zero sometimes become -1E-16, to prevent problems with the logarithm we filter them
    ev = filter(x -> x > 0.0, ev)
    # The von Neumann entropy for the reduced density operator
    entropy = -sum(ev .* log2.(ev))

    return entropy
end


"""
    sample_from_mps!(mps::MPS, gauge_input::Bool=true)::Vector{Int64}

Generate a sample from the probability distribution of basis states defined by the MPS following A. Ferris, G. Vidal, PRB 85 165146 (2021). If the flag gague_input is set to true, the input MPS will be put in right canoncial gauge and normalized, which is a requirement for the algorithm to work. In case one is sure that the MPS is already properly gauged and normalized the gauging step can be spared by setting the flag to false. In this case the input will stay untouched.
"""
function sample_from_mps!(mps::MPS, gauge_input::Bool = true)::Vector{Int64}
    # Extact the length and check if the given position is reasonable
    N = length(mps)
    res = zeros(Int64, N)

    # Put the state into right canonical gauge and make sure it is normalized
    if gauge_input
        gaugeMPS!(mps, :right, true)
    end

    # Now sample from the MPS
    A = 0
    p = 0.0
    M = mps[1][1, :, :]
    for i = 1:N
        d = size(M, 2)
        pacc = 0
        r = rand()
        for l = 1:d
            # Prepare the basis state
            basis_state = zeros(d)
            basis_state[l] = 1.0
            # Contract it into the physical index of the tensor
            A = contract_tensors(M, [2], basis_state, [1])
            # Determine the probability for the basis vector e_l
            p = contract_tensors(conj(A), [1], A, [1])
            pacc += real(p[1])
            if r < pacc
                res[i] = l
                break
            end
        end
        if i < N
            # Contract the result into the next tensor
            M = contract_tensors(1 / sqrt(real(p[1])) * A, [1], mps[i+1], [1])
        end
    end
    return res
end

"""
    sample_from_mps(mps::MPS)::Vector{Int64}

Generate a sample from the probability distribution of basis states defined by the MPS following A. Ferris, G. Vidal, PRB 85 165146 (2021). 
"""
function sample_from_mps(mps::MPS)::Vector{Int64}
    return sample_from_mps!(deepcopy(mps))
end

"""
    svd_compress_mps(mps::MPS, Dmax::Int,, tol::Real = 0.0)::MPS

Compress a given MPS applying an a singular value decomposition at each bond. If Dmax > 0 is supplied (and simultaneously tol = 0.0), a maximum of Dmax singular values is kept, thus truncating the MPS to one with maximum bond dimension Dmax. If tol > 0 is given (and simultaneously Dmax = 0) then all singular values > tol are kept. If both are specified then at most Dmax singular values > tol are kept.
"""
function svd_compress_mps(mps::MPS, Dmax::Int, tol::Real = 0.0)::MPS
    # One of the two parameters has to be larger than zero
    @assert((Dmax > 0) || (tol > 0))
    # Extract the length and prepare a result
    N = length(mps)
    res = deepcopy(mps)
    # Gauge in both directions that redundant dimensions at the boundaries are removed
    gaugeMPS!(res, :left)
    gaugeMPS!(res, :right)
    Dnew = 0
    for i = 1:N-1
        # Check if the given bond dimension is larger than Dmax
        if size(res[i], 2) > Dmax
            Dl1, _, d1 = size(res[i])
            _, Dr2, d2 = size(res[i+1])
            tmp = contract_tensors(res[i], [2], res[i+1], [1])
            tmp = reshape(tmp, (Dl1 * d1, Dr2 * d2))
            U, S, V = svd(tmp)
            M = diagm(S) * V'
            # Now truncate
            if Dmax > 0 && tol == 0.0
                U = U[:, 1:Dmax]
                M = M[1:Dmax, :]
                Dnew = Dmax
            else
                ind = findall(x -> x > tol, S)
                if Dmax > 0 && length(ind) > Dmax
                    ind = ind[1:Dmax]
                end
                U = U[:, ind]
                M = M[ind, :]
                Dnew = length(ind)
            end
            # Reshape and set new tensors
            res[i] = permutedims(reshape(U, (Dl1, d1, Dnew)), (1, 3, 2))
            res[i+1] = reshape(M, (Dnew, Dr2, d2))
        end
    end
    return res
end

"""
    decompose_into_mpo(M::Matrix{T}, d::Vector{Int}) where T <:Number

Given a many-body operator in form of a dense matrix, decompose it into an MPO. It is assumed that the many-body operator follows the index convention of Julia's built-in kronecker product. The vector d specifies the local dimensions of the Hilbert space. 

# Examples
Decomposing a simple two-site operator made up from sum of Pauli terms into MPO form.
```julia-repl
julia>  Id = [1.0 0; 0.0 1.0]
julia>  X = [0.0 1.0; 1.0 0.0]
julia>  Z = [1.0 0.0; 0.0 -1.0]
julia>  H = kron(X,X) + kron(Id,Z) + kron(Z,Id)
julia>  mpo = decompose_into_mpo(H, 2)
```
"""
function decompose_into_mpo(M::Matrix{T}, d::Vector{Int})::MPO{T} where {T<:Number}
    N = length(d)
    dim = prod(d)
    if size(M) != (prod(d), prod(d))
        throw(ArgumentError("matrix not compatible with specified local dimensions, got size(M)=$(repr(size(M))), d=$(repr(d))"))
    end
    # The MPO holding the result
    res = MPO{T}(undef, N)
    # Reshape the matrix and rearrange the indices, such that the row and column index for each site are adjacent and add dummy indices 1 at the boundaries    
    A = reshape(M, (1, reverse(d)..., reverse(d)..., 1))
    ind = collect(Iterators.flatten(zip(collect(N:-1:1), collect(2N:-1:N+1))))
    A = permutedims(A, (1, ind .+ 1..., 2 * length(d) + 2))
    # Now split it into tensors using an SVD
    Dl = 1
    for i = 1:N-1
        # Take the first three indices together
        dims = size(A)
        A = reshape(A, (prod(dims[1:3]), prod(dims[4:end])))
        # SVD the matrix representation
        U, S, V = svd(A)
        A = diagm(S) * V'
        Dr = size(U, 2)
        # Extract the tensor and the remaining part
        U = reshape(U, (Dl, d[i], d[i], Dr))
        res[i] = permutedims(U, (1, 4, 2, 3))
        A = reshape(A, (Dr, dims[4:end]...))
        Dl = Dr
    end
    res[N] = permutedims(A, (1, 4, 2, 3))

    return res
end

"""
    decompose_into_mpo(M::Matrix{T}, d::Int) where T <:Number

Simplified interface to the more general method assuming that all local dimensions are equal to d.
"""
function decompose_into_mpo(M::Matrix{T}, d::Int)::MPO{T} where {T<:Number}
    dl, dr = size(M)
    if dl != dr
        throw(ArgumentError("local dimensions must all be the same, obtained a matrix with dimensions $(repr((dl,dr)))"))
    end
    N = Int(round(log(d, dl)))
    return decompose_into_mpo(M, d * ones(Int64, N))
end