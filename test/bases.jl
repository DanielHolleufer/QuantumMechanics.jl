using QuantumMechanics
using Test

@testset verbose = true "bases" begin
    @testset "GenericBasis" begin
        generic_basis_4 = GenericBasis(4)
        generic_basis_21 = GenericBasis(21)
        generic_basis_5 = GenericBasis(5)
        generic_basis_5f0 = GenericBasis(5.0)

        @test generic_basis_5f0 isa GenericBasis{<:Integer}
        @test generic_basis_5 == generic_basis_5f0
        @test generic_basis_4 == generic_basis_4
        @test generic_basis_21 == generic_basis_21
        @test generic_basis_4 != generic_basis_21

        @test length(generic_basis_4) == 4
        @test length(generic_basis_21) == 21

        @test_throws Exception GenericBasis(0)
        @test_throws Exception GenericBasis(-17)
        @test_throws Exception GenericBasis(5.1)
        @test_throws Exception GenericBasis([2 4])
        @test_throws Exception GenericBasis([2.0 1.1])
        @test_throws Exception GenericBasis([1 2; 3 4])
    end

    @testset "CompositeBasis" begin
        generic_basis_1 = GenericBasis(4)
        generic_basis_2 = GenericBasis(20)
        fock_basis_1 = FockBasis(20)
        fock_basis_2 = FockBasis(50, 25)
        spin_basis_1 = SpinBasis(1//2)
        spin_basis_2 = SpinBasis(2)

        composite_generic = CompositeBasis(generic_basis_1, generic_basis_2)
        composite_fock = CompositeBasis(fock_basis_1, fock_basis_2)
        composite_spin = CompositeBasis(spin_basis_1, spin_basis_2)
        composite_mix_1 = CompositeBasis(generic_basis_1, fock_basis_1, spin_basis_1)
        composite_mix_2 = CompositeBasis(generic_basis_2, fock_basis_2, spin_basis_2)
        composite_mix_all = CompositeBasis(
            generic_basis_1,
            fock_basis_1,
            spin_basis_1,
            generic_basis_2,
            fock_basis_2,
            spin_basis_2,
        )
        composite_nested_1 = CompositeBasis(composite_mix_1, composite_mix_2)
        composite_nested_2 = CompositeBasis(
            (composite_mix_1.bases..., composite_mix_2.bases...),
            [composite_mix_1.dimension; composite_mix_2.dimension],
        )
        composite_mix_generic_1 = CompositeBasis(composite_mix_1, generic_basis_2)
        composite_mix_generic_2 = CompositeBasis(composite_mix_2, generic_basis_1)

        @test composite_generic != generic_basis_1
        @test composite_generic != generic_basis_2
        @test composite_generic != composite_fock
        @test composite_generic != composite_spin
        @test composite_fock != composite_spin
        @test composite_mix_1 != composite_mix_2
        @test composite_mix_all != composite_nested_1
        @test composite_mix_all == composite_nested_2
        @test composite_mix_generic_1 !=
            CompositeBasis(generic_basis_1, fock_basis_1, spin_basis_1, generic_basis_2)
        @test composite_mix_generic_2 !=
            CompositeBasis(generic_basis_2, fock_basis_2, spin_basis_2, generic_basis_1)
        @test composite_mix_generic_1 != composite_mix_generic_2

        @test composite_generic.dimension == [
            generic_basis_1.dimension
            generic_basis_2.dimension
        ]
        @test composite_fock.dimension == [fock_basis_1.dimension; fock_basis_2.dimension]
        @test composite_spin.dimension == [spin_basis_1.dimension; spin_basis_2.dimension]
        @test composite_mix_1.dimension == [
            generic_basis_1.dimension
            fock_basis_1.dimension
            spin_basis_1.dimension
        ]
        @test composite_mix_2.dimension == [
            generic_basis_2.dimension
            fock_basis_2.dimension
            spin_basis_2.dimension
        ]
        @test composite_nested_1.dimension == [
            length(composite_mix_1)
            length(composite_mix_2)
        ]
        @test composite_nested_1.dimension != [
            composite_mix_1.dimension
            composite_mix_2.dimension
        ]
        @test composite_mix_1.dimension != composite_mix_2.dimension
        @test composite_nested_1.dimension != composite_mix_all.dimension

        @test composite_generic.bases == (generic_basis_1, generic_basis_2)
        @test composite_fock.bases == (fock_basis_1, fock_basis_2)
        @test composite_spin.bases == (spin_basis_1, spin_basis_2)
        @test composite_mix_1.bases == (generic_basis_1, fock_basis_1, spin_basis_1)
        @test composite_mix_2.bases == (generic_basis_2, fock_basis_2, spin_basis_2)
        @test composite_mix_1.bases != composite_mix_2.bases

        @test length(composite_generic) == length(generic_basis_1) * length(generic_basis_2)
        @test length(composite_fock) == length(fock_basis_1) * length(fock_basis_2)
        @test length(composite_spin) == length(spin_basis_1) * length(spin_basis_2)
        @test length(composite_mix_1) ==
            length(generic_basis_1) * length(fock_basis_1) * length(spin_basis_1)
        @test length(composite_mix_2) ==
            length(generic_basis_2) * length(fock_basis_2) * length(spin_basis_2)

        @test_throws Exception CompositeBasis(
            (composite_mix_1.bases..., composite_mix_2.bases...),
            [composite_mix_2.dimension; composite_mix_1.dimension],
        )

        tensor_generic = tensor(generic_basis_1, generic_basis_2)
        tensor_fock = tensor(fock_basis_1, fock_basis_2)
        tensor_spin = tensor(spin_basis_1, spin_basis_2)
        tensor_mix_1 = tensor(generic_basis_1, fock_basis_1, spin_basis_1)
        tensor_mix_2 = tensor(generic_basis_2, fock_basis_2, spin_basis_2)
        tensor_mix_12 = tensor(composite_mix_1, composite_mix_2)
        tensor_mix_all = tensor(
            generic_basis_1,
            fock_basis_1,
            spin_basis_1,
            generic_basis_2,
            fock_basis_2,
            spin_basis_2,
        )
        tensor_mix_generic = tensor(composite_mix_1, generic_basis_2)
        tensor_generic_mix = tensor(generic_basis_2, composite_mix_1)

        @test tensor(generic_basis_1) == generic_basis_1
        @test composite_generic == tensor_generic
        @test composite_fock == tensor_fock
        @test composite_spin == tensor_spin
        @test composite_mix_1 == tensor_mix_1
        @test composite_mix_2 == tensor_mix_2
        @test composite_mix_all == tensor_mix_all
        @test tensor_mix_12 == tensor_mix_all
        @test tensor_mix_generic ==
            CompositeBasis(generic_basis_1, fock_basis_1, spin_basis_1, generic_basis_2)
        @test tensor_generic_mix ==
            CompositeBasis(generic_basis_2, generic_basis_1, fock_basis_1, spin_basis_1)
        @test tensor_mix_generic != tensor_generic_mix

        @test tensor_generic == generic_basis_1 ⊗ generic_basis_2
        @test tensor_fock == fock_basis_1 ⊗ fock_basis_2
        @test tensor_spin == spin_basis_1 ⊗ spin_basis_2
        @test tensor_mix_1 == generic_basis_1 ⊗ fock_basis_1 ⊗ spin_basis_1
        @test tensor_mix_1 != generic_basis_1 ⊗ spin_basis_1 ⊗ fock_basis_1
        @test tensor_mix_all ==
            generic_basis_1 ⊗ fock_basis_1 ⊗ spin_basis_1 ⊗ generic_basis_2 ⊗
              fock_basis_2 ⊗ spin_basis_2
    end

    @testset "FockBasis" begin
        fock_basis_10 = FockBasis(10)
        fock_basis_15 = FockBasis(15)
        fock_basis_15_3 = FockBasis(15, 3)
        fock_basis_15f = FockBasis(15.0)
        fock_basis_15_3f = FockBasis(15, 3.0)
        fock_basis_15f_3f = FockBasis(15.0, 3.0)

        @test fock_basis_10.offset == 0

        @test fock_basis_10 == fock_basis_10
        @test fock_basis_15 == fock_basis_15
        @test fock_basis_15 == fock_basis_15f
        @test fock_basis_15_3 == fock_basis_15_3
        @test fock_basis_15_3 == fock_basis_15_3f
        @test fock_basis_15_3 == fock_basis_15f_3f
        @test fock_basis_10 != fock_basis_15
        @test fock_basis_15 != fock_basis_15_3
        @test fock_basis_15 != fock_basis_15_3

        @test length(fock_basis_10) == 11
        @test length(fock_basis_15) == 16
        @test length(fock_basis_15_3) == 13

        @test_throws Exception FockBasis(10, 12)
        @test_throws Exception FockBasis(-3)
        @test_throws Exception FockBasis(10, -3)
        @test_throws Exception FockBasis(-2, -3)
        @test_throws Exception FockBasis(4.1)
        @test_throws Exception FockBasis(4.0, 2.1)
        @test_throws Exception FockBasis(4.1, 2.1)
    end

    @testset "SpinBasis" begin
        spin_basis_1_2 = SpinBasis(1//2)
        spin_basis_1 = SpinBasis(1)
        spin_basis_1_1 = SpinBasis(1//1)
        spin_basis_17_2 = SpinBasis(17//2)

        @test spin_basis_1_2 == spin_basis_1_2
        @test spin_basis_1 == spin_basis_1
        @test spin_basis_1_1 == spin_basis_1_1
        @test spin_basis_1 == spin_basis_1_1
        @test spin_basis_1_2 != spin_basis_1
        @test spin_basis_1_2 != spin_basis_17_2

        @test length(spin_basis_1_2) == 2
        @test length(spin_basis_1) == 3
        @test length(spin_basis_17_2) == 18

        @test_throws Exception SpinBasis(-1)
        @test_throws Exception SpinBasis(-1//2)
        @test_throws Exception SpinBasis(-1//1)
        @test_throws Exception SpinBasis(0)
        @test_throws Exception SpinBasis(0//2)
        @test_throws Exception SpinBasis(0//0)
        @test_throws Exception SpinBasis(1//3)
        @test_throws Exception SpinBasis(0.25)
    end
end
