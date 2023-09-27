using QuantumMechanics
using Test

@testset verbose = true "StateVectors" begin
    @testset "Bra and Ket" begin
        u = [0.0 + 0.0im, 0.0 - 0.0im, 0.0 + 1.0im / 2, 0.0 + 0.0im]
        v = [1.0 / sqrt(2.0) + 0.0im, 0.0 - 1.0im / 2.0, 0.0 + 0.0im, -1.0 / 2.0 + 0.0im]
        basis_1 = GenericBasis(4)
        bra_u = Bra(basis_1, u)
        ket_u = Ket(basis_1, u)
        bra_v = Bra(basis_1, v)
        ket_v = Ket(basis_1, v)
        @test bra_u == bra_u
        @test ket_u == ket_u
        @test bra_u != ket_u
        @test bra_u != bra_v
        @test ket_u != ket_v
        @test bra_u != ket_v
        @test ket_u != bra_v
        @test bra_u.data == ket_u.data'
        @test bra_u.data != ket_u.data

        bra_u_ad = Bra(basis_1, u')
        ket_u_ad = Ket(basis_1, u')
        bra_u_conj = Bra(basis_1, conj(u))
        ket_v_conj = Ket(basis_1, conj(u))
        @test bra_u_ad == bra_u
        @test ket_u_ad == ket_u
        @test bra_u_conj != bra_u
        @test ket_v_conj != ket_u

        @test Bra{typeof(basis_1)}(basis_1, u) == bra_u
        @test Ket{typeof(basis_1)}(basis_1, u) == ket_u
        @test Bra{typeof(basis_1)}(basis_1, u') == bra_u
        @test Ket{typeof(basis_1)}(basis_1, u') == ket_u

        basis_2 = GenericBasis(5)
        @test_throws Exception Bra(basis_2, u)
        @test_throws Exception Ket(basis_2, u)
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
        @test bra_u * ket_u ≈ 1.0 atol = 1.0e-12
        @test bra_v * ket_v ≈ 1.0 atol = 1.0e-12
        @test bra_u * ket_v ≈ -im * sqrt(3 / 2) / 2 atol = 1.0e-12
        @test bra_v * ket_u ≈ im * sqrt(3 / 2) / 2 atol = 1.0e-12

        @test_throws Exception bra_u + ket_v
        @test_throws Exception ket_u + bra_v
        @test_throws Exception bra_u - ket_v
        @test_throws Exception ket_u - bra_v
        @test_throws Exception b / bra_u
        @test_throws Exception b / ket_u

        w = [0 + im / sqrt(3), 1 / sqrt(3) + 0im, -1 / sqrt(3) + 0im, 0 + 0im]
        basis_2 = FockBasis(3)
        bra_w = Bra(basis_2, w)
        ket_w = Ket(basis_2, w)
        @test_throws Exception ket_u + ket_w
        @test_throws Exception bra_u + bra_w
        @test_throws Exception ket_u - ket_w
        @test_throws Exception bra_u - bra_w
        @test_throws Exception bra_u * ket_w
        @test_throws Exception bra_w * ket_u
    end

    @testset "copy" begin
        v = [1.0 / sqrt(2.0) + 0.0im, 0.0 - 1.0im / 2.0, 0.0 + 0.0im, -1.0 / 2.0 + 0.0im]
        basis = GenericBasis(4)
        bra_v = Bra(basis, v)
        ket_v = Ket(basis, v)
        bra_v_copy = copy(bra_v)
        ket_v_copy = copy(ket_v)
        @test bra_v == bra_v_copy
        @test ket_v == ket_v_copy

        bra_v_copy.data[1] = -1.0im / sqrt(2.0)
        ket_v_copy.data[1] = 1.0im / sqrt(2.0)
        @test bra_v != bra_v_copy
        @test ket_v != ket_v_copy

        v_copy = copy(v)
        v_copy[1] = 1.0im / sqrt(2.0)
        @test bra_v_copy == Bra(basis, v_copy)
        @test ket_v_copy == Ket(basis, v_copy)
    end

    @testset "dagger" begin
        v_1 = [1.0 / sqrt(2.0) + 0.0im, 0.0 - 1.0im / 2.0, 0.0 + 0.0im, -1.0 / 2.0 + 0.0im]
        basis_1 = GenericBasis(4)
        @test dagger(Ket(basis_1, v_1)) == Bra(basis_1, v_1)
        @test dagger(Bra(basis_1, v_1)) == Ket(basis_1, v_1)
        @test dagger(Ket(basis_1, v_1)) == Bra(basis_1, v_1')
        @test dagger(Bra(basis_1, v_1)) == Ket(basis_1, v_1')
        @test dagger(dagger(Ket(basis_1, v_1))) == Ket(basis_1, v_1)
        @test dagger(dagger(Bra(basis_1, v_1))) == Bra(basis_1, v_1)
        @test dagger(Ket(basis_1, v_1)) != Bra(basis_1, conj(v_1))
        @test dagger(Bra(basis_1, v_1)) != Ket(basis_1, conj(v_1))

        @test Ket(basis_1, v_1)' == adjoint(Ket(basis_1, v_1))
        @test Bra(basis_1, v_1)' == adjoint(Bra(basis_1, v_1))
        @test Ket(basis_1, v_1)' == Bra(basis_1, v_1')
        @test Bra(basis_1, v_1)' == Ket(basis_1, v_1')
        @test Ket(basis_1, v_1)'' == Ket(basis_1, v_1)
        @test Bra(basis_1, v_1)'' == Bra(basis_1, v_1)
        @test Ket(basis_1, v_1)' != Bra(basis_1, conj(v_1))
        @test Bra(basis_1, v_1)' != Ket(basis_1, conj(v_1))
    end

    @testset "norm and normalize" begin
        v = [1.0 / sqrt(2.0) + 0.0im, 0.0 - 1.0im / 2.0, 0.0 + 0.0im, -1.0 / 2.0 + 0.0im]
        w = [1.0 + 0.0im, 0.0 - 1.0im / 2.0, 0.0 + 2.0im, -1.0 / 2.0 + 1.0im / 2.0]
        basis = GenericBasis(4)
        bra_v = Bra(basis, v)
        ket_v = Ket(basis, v)
        bra_w = Bra(basis, w)
        ket_w = Ket(basis, w)

        @test norm(bra_v) ≈ 1.0 atol = 1.0e-12
        @test norm(ket_v) ≈ 1.0 atol = 1.0e-12
        @test norm(bra_w) ≈ sqrt(1.0 + 0.25 + 4.0 + 0.5) atol = 1.0e-12
        @test norm(ket_w) ≈ sqrt(1.0 + 0.25 + 4.0 + 0.5) atol = 1.0e-12
        @test norm(normalize(bra_w)) ≈ 1.0 atol = 1.0e-12
        @test norm(normalize(ket_w)) ≈ 1.0 atol = 1.0e-12
        
        normalize!(bra_w)
        normalize!(ket_w)
        @test norm(bra_w) ≈ 1.0 atol = 1.0e-12
        @test norm(ket_w) ≈ 1.0 atol = 1.0e-12
    end
end
