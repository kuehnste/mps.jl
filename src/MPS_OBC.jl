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
function random_mps_obc(N::Int, D::Int, d, tensortype::Type{T}=ComplexF64)::MPS{T} where T
    # Ensure the the system size is at least two
    @assert(N>1)

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

    # Left tensor (row vector)
    mps[1] = rand(tensortype, 1, D, dim[1])

    # Right tensor (column vector)
    mps[N] = rand(tensortype, D, 1, dim[N])

    # Tensors in between
    for i = 2:N - 1
        mps[i] = rand(tensortype, D, D, dim[i])
    end

    return mps
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
        overlap = contract_tensors(conj(mps1[i]), [1;3], overlap, [1;3])
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
        val = contract_tensors(val, [1;4], mpo[i], [1;3])
        val = contract_tensors(val, [1;4], mps[i], [1;3])
    end

    return val[1]    
end

"""
    gaugeMPS(mps::MPS{T}, direction::Symbol=:right, normalize::Bool=false)::MPS{T}
    
Function to bring an MPS with OBC in canonical form
direction: tells, whether it will be left or right canonical gauge
If normalize is set to true, the resulting state will be normalized.
"""
function gaugeMPS(mps::MPS{T}, direction::Symbol=:right, normalize::Bool=false)::MPS{T} where T    
    mps_gauged = deepcopy(mps)
    gaugeMPS!(mps, direction, normalize)    
    return mps_gauged
end

"""
    gaugeMPS!(mps::MPS,direction::Dir=right,)
    
Bring an MPS with OBC in canonical form. Direction tells, whether it will be left or right canonical gauge If normalize is set to true, the resulting state will be
normalized. This function overwrites the input MPS with its gauged version
"""
function gaugeMPS!(mps::MPS, direction::Symbol=:right, normalize::Bool=false)
    # Check that we got a meaningful direction
    if (direction!=:left && direction!=right)
        throw(ArgumentError("direction must be :left or :right, got $(repr(direction))"))
    end
    # Start gauging
    N = length(mps)
    if (direction == :left)
        # Case that I want left canonical gauge
        M = mps[1]
        for i = 1:N - 1
            mps[i], res = gauge_site(M, direction)
            M = contract_tensors(res, [2], mps[i + 1], [1])
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
            M = contract_tensors(mps[i - 1], [2], res, [1], [1;3;2])
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
function gauge_site(A::Site{T}, direction::Symbol)::Tuple{Site{T},Matrix{T}} where T
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
    find_groundstate(H::Array{Operator},D::Int64,acc::Float64,max_sweeps::Int64=-1)

Given a Hamiltonian MPO H, find an MPS approximation for its ground state with bond 
dimension D converged to a relative accuracy acc. If a positive number for max_sweeps
is given, the maximum amount of iterations is limited to max_sweeps
"""
function find_groundstate(H::Vector{Operator{T}}, D::Int64, d::Int64, acc::Float64, max_sweeps::Int64) where T
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
        for i = 1:N - 1
	        Left = LR[i]
            Right = LR[i + 1]
            E_local, M = solve_eigenvalue_problem(Left, Right, H[i])
            mps[i], _ = gauge_site(M, :left)        
            LR[i + 1] = update_left(Left, mps[i], mps[i], H[i])
        end 

        # From left to right starting by 1 up to N-1
        for i = N:-1:2
            Left = LR[i]
            Right = LR[i + 1]
            E_local, M = solve_eigenvalue_problem(Left, Right, H[i])
            mps[i], res = gauge_site(M, :right)
            LR[i] = update_right(Right, mps[i], mps[i], H[i]) 
        end
      
      # Check convergence criterion
        if (abs((Eold - E_local) / Eold) < acc)    
	        E0 = E_local
	        # Restore normalisation (permutation is needed to bring the indices back into right order after contracting)
	        mps[1] = contract_tensors(mps[1], [2], res, [1], [1;3;2])  
	        println("Convergence to desired accuracy achieved")
	        break
        elseif (max_sweeps > 0 && num_of_sweeps >= max_sweeps)
	        E0 = E_local
	        # Restore normalisation
	        mps[1] = contract_tensors(mps[1], [2], res, [1], [1;3;2]) 
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
    Htemp = permutedims(Htemp, [1;5;3;2;6;4])
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
    Htemp = permutedims(Htemp, [1;5;3;2;6;4])
    dim1, dim2, dim3, dim4, dim5, dim6 = size(Htemp)
    Heff = reshape(Htemp, (dim1 * dim2 * dim3, dim4 * dim5 * dim6))

    E_local, M = eigs(Heff, nev=1, which=:SR)
    M = reshape(M, (dim1, dim5, dim3))

    return real(E_local[1]), M
end

"""
    setup_R(H::MPO{T1},mps::MPS{T2}) where {T1,T2}
    
Function to calculate the partial contractions needed to form the effective Hamiltonian
"""
function setup_R(H::MPO{T1}, mps::MPS{T2}) where {T1,T2}
    N = length(H)
    if (T1 <: T2)
        entrytype = T2
    else
        entrytype = T1
    end
    LR = Vector{Array{entrytype,3}}(undef, N + 1)

    # We need only N-1 partial contractions, however, we set the edges to dummy values 1 that we can recoursivly compute every computation resuing the others
    LR[1] = ones(entrytype, 1, 1, 1)
    LR[N + 1] = ones(entrytype, 1, 1, 1)

    # Now compute the partial contractions starting from the right (as I start sweeping on the left in the ground state search)
    for i = N:-1:2
        tmp = contract_tensors(LR[i + 1], [3], mps[i], [2])
        tmp = contract_tensors(H[i], [2;4], tmp, [2;4])
        LR[i] = contract_tensors(conj(mps[i]), [2;3], tmp, [3;2])
    end

    return LR
end

"""
  update_left(LR, M_left::Site, M_right::Site, W::Operator)
    
Function to perform an update the partial contractions starting from the left stored in LR 
"""
function update_left(LR, M_left::Site, M_right::Site, W::Operator)
    res = contract_tensors(LR, [3], M_right, [1])
    res = contract_tensors(W, [1;4], res, [2;4])
    res = contract_tensors(conj(M_left), [1;3], res, [3;2])

    return res
end

"""
  update_right(LR, M_left::Site, M_right::Site, W::Operator)
    
Function to perform an update the partial contractions starting from the left stored in LR 
"""
function update_right(LR, M_left::Site, M_right::Site, W::Operator)
    res = contract_tensors(LR, [3], M_right, [2])
    res = contract_tensors(W, [2;4], res, [2;4])
    res = contract_tensors(conj(M_left), [2;3], res, [3;2])  

    return res
end

"""
  apply_operator(operator::MPO{T1}, mps::MPS{T2})::MPS where {T1,T2}
    
Apply an opterator given as MPO to an MPS
"""
function apply_operator(operator::MPO{T1}, mps::MPS{T2})::MPS where {T1,T2}
    N = length(mps)
  
    # Generate a new MPS of the correct type
    if (T1 <: T2)
        res = MPS{T2}(undef, N)
    else
        res = MPS{T1}(undef, N)
    end     

    # Apply Hamilton to MPS and generate new MPS
    for i = 1:N
        temp = contract(operator[i], [4], mps[i], [3])
        temp = permute(temp, [1 4 2 5 3])
        dim1, dim2, dim3, dim4, dim5 = size(temp)
        res[i] = reshape(temp, [dim1 * dim2, dim3 * dim4, dim5])
    end
    return res
end

"""
   approximate_mps(D::Int64,mps0::MPS{T},acc::Float64)::MPS{T}

Function to calculate an approximate MPS with bond dimension D to a given MPS mps0 such that ||mps - mps0||_2 -> min
"""
function approximate_mps(D::Int64, mps0::MPS{T}, acc::Float64)::MPS{T} where T
    # Number of spins
    N = length(mps0)
    # Physical dimension
    d = size(mps0{1}, 3)

    # Initial guess for MPS
    mps = random_mps_obc(N, D, d, T)
    # Right normlised form
    gaugeMPS!(mps, :right)

    # Initialise L and R
    LR = initialise(mps0, mps)

    # Optimization
    num_of_sweeps = 0
    dev_old = 0
    alldevs = zeros(2 * N - 2, 1)

    while (1)    
        num_of_sweeps = num_of_sweeps + 1
        println("Sweep ", num_of_sweeps)
        count = 1

        # Sweep form left to right
        for i = 1:N - 1
	      Left = LR[i]
	      Right = LR[i + 1]
	      temp = contract_tensors(Right, [1], mps0[i], [2])
	      mps[i] = contract_tensors(Left, [2], temp, [2])
	      alldevs[count] = check_accu(mps[i])
	      _, mps[i] = gauge_site(mps[i], :left)  
	      LR[i + 1] = update_left(Left, mps[i], mps0[i])
	      count = count + 1
        end

        # Sweep form right to left
        for i = N:-1:2
            Left = LR[i]
            Right = LR[i + 1]
            temp = contract_tensors(Right, [1], mps0[i], [2])
            mps[i] = contract_tensors(Left, [2], temp, [2])
            alldevs[count] = check_accu(mps[i])
            A, mps[i] = gauge_site(mps[i], :right)    
            LR[i] = update_right(Right, mps0[i], mps[i])
            count = count + 1
        end

        if (std(alldevs) / abs(mean(alldevs)) < acc)
	        break
        end      
    end
    # Restore normalisation
    mps[1] = permute(contract_tensors(mps[1], [2], A, [1]), [1;3;2])
    
    return mps
end

