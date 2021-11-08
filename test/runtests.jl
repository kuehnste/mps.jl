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
