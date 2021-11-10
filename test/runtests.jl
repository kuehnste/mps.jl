using Test 
using LinearAlgebra
using MatrixProductStates

##############################################################
#   Test the generation and normalization of an random MPS   #
##############################################################

@testset "Random MPS" begin
    # Use left canonical form
    for N = 10:2:20
        mps_real = random_mps_obc(N, 5, 3, Float64)
        mps_complex = random_mps_obc(N, 5, 3, ComplexF64)
        gaugeMPS!(mps_real, :left, true)
        gaugeMPS!(mps_complex, :left, true)
        @test isapprox(calculate_overlap(mps_real, mps_real), 1.0)
        @test isapprox(calculate_overlap(mps_complex, mps_complex), 1.0 + 0.0im)
    end
    # Use right canonical form
    for N = 10:2:20
        mps_real = random_mps_obc(N, 5, 3, Float64)
        mps_complex = random_mps_obc(N, 5, 3, ComplexF64)
        gaugeMPS!(mps_real, :right, true)
        gaugeMPS!(mps_complex, :right, true)
        @test isapprox(calculate_overlap(mps_real, mps_real), 1.0)
        @test isapprox(calculate_overlap(mps_complex, mps_complex), 1.0 + 0.0im)
    end
end

@testset "Identity operator" begin    
    for N = 10:2:20
        for d = 2:5
            mps_complex = random_mps_obc(N, 5, d, ComplexF64)
            gaugeMPS!(mps_complex, :left, true)
            id_mpo = getIdentityMPO(N, d)
            mps_operator_applied = apply_operator(id_mpo, mps_complex)
            @test isapprox(calculate_overlap(mps_complex, mps_operator_applied), 1.0 + 0.0im)  
        end      
    end    
end

@testset "Canonical form" begin
    mps = random_mps_obc(10, 5, 2, ComplexF64)
    mps_left_gauged = gaugeMPS(mps, :left, true)
            mps_right_gauged = gaugeMPS(mps, :right, true)
    # Check the left canonical gauge
    for i = 1:length(mps_left_gauged)
        Dr = size(mps_left_gauged[i], 2)
        res = contract_tensors(conj(mps_left_gauged[i]), [1;3], mps_left_gauged[i], [1;3])
        @test isapprox(res, Matrix((1.0 + 0.0im) * I, Dr, Dr))
    end
    # Check the right canonical gauge
    for i = 1:length(mps_right_gauged)
        Dl = size(mps_right_gauged[i], 1)
        res = contract_tensors(conj(mps_right_gauged[i]), [2;3], mps_right_gauged[i], [2;3])
        @test isapprox(res, Matrix((1.0 + 0.0im) * I, Dl, Dl))
    end
end

@testset "Inplace operator application" begin
    # Check that we get an error, if we try to overwrite a real MPS with a complex MPO
    mps = random_mps_obc(10, 9, 2, Float64)
    H = getHeisenbergMPO(10, 1, 1)
    @test_throws InexactError apply_operator!(H, mps)

    # Compute result of the contraction and store in a new MPS, compare it with result of inplace contraction
    for N = 10:2:20
        # Prepare a random MPS and gauge it
        mps = random_mps_obc(N, 15, 2, ComplexF64)
        gaugeMPS!(mps, :left, true)
        # Some MPOs for testing
        H = getIsingMPO(N, 1.0, 0.9)
        # Apply the Hamiltonian generating a copy
        mps_new = apply_operator(H, mps)
        gaugeMPS!(mps_new, :left, true)
        # Now in place
        apply_operator!(H, mps)
        gaugeMPS!(mps, :left, true)
        # Check the expectation values
        @test isapprox(calculate_overlap(mps, mps_new), 1.0 + 0.0im)
    end
end

@testset "Summing MPOs" begin
    for N = 10:2:20
        # Prepare a random MPS and gauge it
        mps = random_mps_obc(N, 15, 2, ComplexF64)
        gaugeMPS!(mps, :left, true)
        # Some MPOs for testing
        id_mpo = getIdentityMPO(N, 2)
        H = getIsingMPO(N, 1.0, 0.9)
        # Compute sums thereof
        O1 = sum_operators(id_mpo, H)
        O2 = sum_operators(H, H)
        # Check the expectation values
        @test isapprox(expectation_value(mps, O1), 1.0 + expectation_value(mps, H))
        @test isapprox(expectation_value(mps, O2), 2.0 * expectation_value(mps, H))
    end
end

@testset "Test entropy computation" begin
    # A product state should always have vanishing entropy
    mps = random_mps_obc(10, 1, 2)
    gaugeMPS!(mps, :left, true)
    for i = 1:length(mps)        
        @test abs(compute_entropy(mps, i)) < 1E-12
    end
    @test abs(compute_entropy(mps)) < 1E-12

    # Now prepare the GHZ state
    all_zeros = product_state_obc(BitArray(zeros(Int8, 10)))
    all_ones = product_state_obc(BitArray(ones(Int8, 10)))
    mps = sum_states(all_zeros, all_ones)
    gaugeMPS!(mps, :left, true)
    for i = 1:length(mps) - 1
        @test isapprox(compute_entropy(mps, i), 1.0)        
    end

    # Compare the entropy computation with the default argument to what one would expect
    mps = random_mps_obc(11, 4, 3)
    entropy1 = compute_entropy(mps)
    entropy2 = compute_entropy(mps, 6)
    @test isapprox(entropy1, entropy2)

    mps = random_mps_obc(8, 6, 3)
    entropy1 = compute_entropy(mps)
    entropy2 = compute_entropy(mps, 4)
    @test isapprox(entropy1, entropy2)
end
