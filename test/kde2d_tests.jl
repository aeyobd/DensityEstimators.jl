
@testset "test kernel2d_epanechnikov" begin
    @test Arya.kernel2d_epanechnikov(0, 0) ≈ 2/π
    @test Arya.kernel2d_epanechnikov(1, 0) ≈ 0
    @test Arya.kernel2d_epanechnikov(-0.7, 0.8) ≈ 0

end


@testset "kernel2d_epanechnikon normalization" begin
    N = 1000000
    x = -1 .+ 2rand(N)
    y = -1 .+ 2rand(N)
    @test sum(Arya.kernel2d_epanechnikov.(x, y)) ≈ N / 4 rtol=1e-2

end


@testset "kernel2d_epanechnikon integration" begin
    @test Arya.integrate2d(Arya.kernel2d_epanechnikov, -1, 1, -1, 1) ≈ 1 rtol=1e-2 broken = true

end
