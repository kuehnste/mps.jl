###########################################################################
#       Basic functionality for contracting and manipulating tensors      #
###########################################################################

using LinearAlgebra


"""
    contract_tensors(T1,cT1::Vector{Int64},T2,cT2::Vector{Int64},ordering::Vector{Int}=Vector{Int64}())

Function to contract to tensors T1 and T2, the contraction is performed as matrix multiplication by reordering the indices and reshaping the tensors
num_of_indices_T1: number of indices of tensor T1
num_of_indices_T2: number of indices of tensor T2
iT1: vector of indices of T1 which are contracted
iT2: vector of indices of T2 which are contracted
"""
function contract_tensors(T1, cT1::Vector{Int}, T2, cT2::Vector{Int}, ordering::Vector{Int}=Vector{Int64}())
    # Make sure the number of indices to contract is the same
    @assert(length(cT1) == length(cT2))

    # Determine size of the tensors
    dimT1 = size(T1)
    dimT2 = size(T2)

    # Make sure the dimensions along the contracted indices are equal
    if (dimT1[cT1] != dimT2[cT2])
        error("Dimensions along the contracted indices must be the same, failed to contract the indices")
    end

    # Get the uncontracted indices of T1 and T2
    uT1 = setdiff(collect(1:ndims(T1)), cT1)
    uT2 = setdiff(collect(1:ndims(T2)), cT2)

    # Bring the contracted indices of T1 to the end
    T1 = permutedims(T1, (uT1..., cT1...))
    # Bring the contracted indices of T2 to the front
    T2 = permutedims(T2, (cT2..., uT2...))

    # Build single index out of all contracted and uncontracted indices to
    # obtain a matrix with two indices
    il = prod(dimT1[uT1])
    ir = prod(dimT2[uT2])
    T1 = reshape(T1, (il, prod(dimT1[cT1])))
    T2 = reshape(T2, (prod(dimT2[cT2]), ir))
    # Actual contraction as matrix multiplication
    res = T1 * T2
  
    # Expand the single indices
    res = reshape(res, (dimT1[uT1]..., dimT2[uT2]...))

    # Permute order if necessary
    if (!isempty(ordering))
        res = permutedims(res, (ordering...,))
    end

    return res
end


