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

        u_2 = [1 / sqrt(2) -1 / sqrt(2)]
        v_2 = [1 0]
        w_2 = [1.0 / sqrt(2), 1.0im / sqrt(2)]
        basis_2 = GenericBasis(2)
        ϕ_2 = Bra(basis_2, u_2)
        ψ_2 = Ket(basis_2, v_2)
        @test ψ_2 == ψ_2
        @test ϕ_2 == ϕ_2
        @test ψ_2 != ϕ_2
        @test typeof(Bra(basis_2, u_2)) == typeof(Bra(basis_2, v_2))
        @test typeof(Bra(basis_2, u_2)) == typeof(Bra(basis_2, w_2))
        @test typeof(Ket(basis_2, u_2)) == typeof(Ket(basis_2, v_2))
        @test typeof(Ket(basis_2, u_2)) == typeof(Ket(basis_2, w_2))
        @test typeof(Bra(basis_2, u_2)) != typeof(Ket(basis_2, w_2))

        @test_throws Exception Bra(basis_1, u_2)
        @test_throws Exception Bra(basis_2, u_1)
        @test_throws Exception Ket(basis_1, u_2)
        @test_throws Exception Ket(basis_2, u_1)
        @test_throws Exception Bra(
            basis_1, (0.0 + 0.0im, 0.0 - 0.0im, 0.0 + 1.0im / 2, 0.0 + 0.0im)
        )
    end

    @testset "Arithmetic and vector operations" begin
        a = exp(-im * π / sqrt(2))
        b = sqrt(2)
        u = [0.0 + im * sqrt(3.0) / 2.0, 0.0 + 0.0im, 0.0 + 1.0im / 2, 0.0 + 0.0im]
        v = [1.0 / sqrt(2.0) + 0.0im, 0.0 - 1.0im / 2.0, 0.0 + 0.0im, -1.0 / 2.0 + 0.0im]
        basis = GenericBasis(4)
        bra_u = Bra(basis, u)
        ket_u = Ket(basis, u)
        bra_v = Bra(basis, v)
        ket_v = Ket(basis, v)

        @test bra_u + bra_v == Bra(basis, u + v)
        @test ket_u + ket_v == Ket(basis, u + v)
        @test bra_u - bra_v == Bra(basis, u - v)
        @test ket_u - ket_v == Ket(basis, u - v)
        @test a * bra_u == Bra(basis, a * u)
        @test a * ket_u == Ket(basis, a * u)
        @test bra_u * a == Bra(basis, u * a)
        @test ket_u * a == Ket(basis, u * a)
        @test bra_u / b == Bra(basis, u / b)
        @test ket_u / b == Ket(basis, u / b)
        @test bra_u * ket_u ≈ 1-0 atol = 1.0e-12
        @test bra_v * ket_v ≈ 1-0 atol = 1.0e-12
        @test bra_u * ket_v ≈ -im * sqrt(3 / 2) / 2 atol = 1.0e-12
        @test bra_v * ket_u ≈ im * sqrt(3 / 2) / 2 atol = 1.0e-12

        @test_throws Exception bra_u + ket_v
        @test_throws Exception ket_u + bra_v
        @test_throws Exception bra_u - ket_v
        @test_throws Exception ket_u - bra_v
        @test_throws Exception b / bra_u
        @test_throws Exception b / ket_u

        w = [0.0 + im / sqrt(3.0), 1.0 / sqrt(3.0) + 0.0im, -1.0 / sqrt(3.0) + 0.0im]
        basis_2 = GenericBasis(3)
        bra_w = Bra(basis_2, w)
        ket_w = Ket(basis_2, w)
        @test_throws Exception ket_u + ket_w
        @test_throws Exception bra_u + bra_w
        @test_throws Exception ket_u - ket_w
        @test_throws Exception bra_u - bra_w
        @test_throws Exception bra_u * ket_w
        @test_throws Exception bra_w * ket_u
    end

    @testset "dagger" begin
        v_1 = [1.0 / sqrt(2.0) + 0.0im, 0.0 - 1.0im / 2.0, 0.0 + 0.0im, -1.0 / 2.0 + 0.0im]
        basis_1 = GenericBasis(4)
        @test dagger(Ket(basis_1, v_1)) == Bra(basis_1, v_1')
        @test dagger(Bra(basis_1, v_1)) == Ket(basis_1, v_1')
        @test dagger(dagger(Ket(basis_1, v_1))) == Ket(basis_1, v_1)
        @test dagger(dagger(Bra(basis_1, v_1))) == Bra(basis_1, v_1)
        @test dagger(Ket(basis_1, v_1)) != Bra(basis_1, v_1)

        @test Ket(basis_1, v_1)' == adjoint(Ket(basis_1, v_1))
        @test Bra(basis_1, v_1)' == adjoint(Bra(basis_1, v_1))
        @test Ket(basis_1, v_1)' == Bra(basis_1, v_1')
        @test Bra(basis_1, v_1)' == Ket(basis_1, v_1')
        @test Ket(basis_1, v_1)'' == Ket(basis_1, v_1)
        @test Bra(basis_1, v_1)'' == Bra(basis_1, v_1)
        @test Ket(basis_1, v_1)' != Bra(basis_1, v_1)
    end
end
