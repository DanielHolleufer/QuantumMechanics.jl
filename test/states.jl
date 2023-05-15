using QuantumMechanics
using Test

@testset "StateVectors" begin
    v = [1.0 / sqrt(2.0) + 0.0im, 0.0 - 1.0im / 2.0, 0.0 + 0.0im, -1.0 / 2.0 + 0.0im]
    b = GenericBasis(4)
    ψ = Ket(b, v)
    @test ψ == ψ

    ϕ = Bra(b, v)
    @test ϕ == ϕ
    @test ψ != ϕ

    @test ϕ * ψ ≈ 1.0 atol = 1.0e-12


    u = [0.0 + 1.0im / sqrt(3.0) 0.0 - 0.0im exp(-2.0im) * sqrt(2 / 3)]
    b = GenericBasis(3)
    ψu = Ket(b, v)
    @test ψu == ψu

    ϕu = Bra(b, v)
    @test ϕu == ϕu
    @test ψu != ϕu

    @test ϕu * ψu ≈ 1.0 atol = 1.0e-12
end
