using QuantumMechanics
using Test

@testset verbose = true "bases" begin
    @testset "GenericBasis" begin
        generic_basis_4 = GenericBasis(4)
        generic_basis_21 = GenericBasis(21)

        @test generic_basis_4 == generic_basis_4
        @test generic_basis_21 == generic_basis_21
        @test generic_basis_4 != generic_basis_21

        @test length(generic_basis_4) == 4
        @test length(generic_basis_21) == 21

        @test_throws Exception GenericBasis(0)
        @test_throws Exception GenericBasis(-17)
    end

    @testset "FockBasis" begin
        fock_basis_10 = FockBasis(10)
        fock_basis_15 = FockBasis(15)
        fock_basis_15_3 = FockBasis(15, 3)

        @test fock_basis_10.offset == 0

        @test fock_basis_10 == fock_basis_10
        @test fock_basis_15 == fock_basis_15
        @test fock_basis_15_3 == fock_basis_15_3
        @test fock_basis_10 != fock_basis_15
        @test fock_basis_15 != fock_basis_15_3
        @test fock_basis_15 != fock_basis_15_3

        @test length(fock_basis_10) == 11
        @test length(fock_basis_15) == 16
        @test length(fock_basis_15_3) == 13

        @test_throws Exception FockBasis(10, 12)
        @test_throws Exception GenericBasis(-3)
    end
end
