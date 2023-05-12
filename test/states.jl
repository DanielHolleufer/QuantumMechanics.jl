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
end