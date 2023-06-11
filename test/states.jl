using QuantumMechanics
using Test

@testset verbose = true "StateVectors" begin
    @testset "Bra and Ket" begin
        u_1 = [0.0 + 0.0im, 0.0 - 0.0im, 0.0 + 1.0im / 2, 0.0 + 0.0im]
        v_1 = [1.0 / sqrt(2.0) + 0.0im, 0.0 - 1.0im / 2.0, 0.0 + 0.0im, -1.0 / 2.0 + 0.0im]
        basis_1 = GenericBasis(4)
        ϕ_1 = Bra(basis_1, u_1)
        ψ_1 = Ket(basis_1, v_1)
        @test ψ_1 == ψ_1
        @test ϕ_1 == ϕ_1
        @test ψ_1 != ϕ_1
        @test ϕ_1 * ψ_1 ≈ 0.0 atol = 1.0e-12

        u_2 = [1 / sqrt(2) -1 / sqrt(2)]
        v_2 = [1 0]
        w_2 = [1.0 / sqrt(2), 1.0im / sqrt(2)]
        basis_2 = GenericBasis(2)
        ϕ_2 = Bra(basis_2, u_2)
        ψ_2 = Ket(basis_2, v_2)
        @test ψ_2 == ψ_2
        @test ϕ_2 == ϕ_2
        @test ψ_2 != ϕ_2
        @test ϕ_2 * ψ_2 ≈ 1 / sqrt(2) atol = 1.0e-12
        @test typeof(Bra(basis_2, u_2)) == typeof(Bra(basis_2, v_2))
        @test typeof(Bra(basis_2, u_2)) == typeof(Bra(basis_2, w_2))
        @test typeof(Ket(basis_2, u_2)) == typeof(Ket(basis_2, v_2))
        @test typeof(Ket(basis_2, u_2)) == typeof(Ket(basis_2, w_2))
        @test typeof(Bra(basis_2, u_2)) != typeof(Ket(basis_2, w_2))

        @test_throws Exception Bra(basis_1, u_2)
        @test_throws Exception Bra(basis_2, u_1)
        @test_throws Exception Bra(
            basis_1, (0.0 + 0.0im, 0.0 - 0.0im, 0.0 + 1.0im / 2, 0.0 + 0.0im)
        )
    end
end
