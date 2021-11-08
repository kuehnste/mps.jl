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
