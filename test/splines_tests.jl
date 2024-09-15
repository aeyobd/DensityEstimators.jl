import DensityEstimators as DE

@testset "spline construction" begin

    t = [0.0, 1.0, 2.0, 3.0, 4.0]
    α = [2, 3]
    k = 3
    s = DE.BSpline(t, α, k)

    @test DE.knots_of(s) == t
    @test DE.coefficients_of(s) == α
    @test DE.order_of(s) == k


    @test_throws DimensionMismatch DE.BSpline(t, α, 0)
    @test_throws DimensionMismatch DE.BSpline(t, α[2:end], k)
    @test_throws DimensionMismatch DE.BSpline(t[2:end], α, k)


    t = [0.0, 1.0, 2.0, 3.0, 4.0]
    α = [2, 3, 4]
    k = 2
    s = DE.BSpline(t, α, k)

    @test_throws ArgumentError DE.BSpline(t[end:-1:1], α, k)
end



@testset "evaluation" begin

    @testset "order 1" begin
        # order 1 spline is just a piecewise constant function

        t = [-0.5, 0.6, 0.8, 1.52]
        α = [    2,    π, -0.3]
        k = 1

        s = DE.BSpline(t, α, k)

        @test s(-1.) ≈ 0
        @test s([1.53, 28]) ≈ [0, 0]
        @test s(0.) ≈ 2
        @test s(0.8) ≈ -0.3
        @test s(-0.5) ≈ 2
        @test s(0.6) ≈ π
        @test s(1.52) ≈ -0.3
        @test s(0.7) ≈ π
        @test s(LinRange(-0.5, 0.59, 100)) ≈ fill(2, 100)
    end

    @testset "order 2" begin

    end
end



@testset "integration" begin

    @testset "quadgk" begin
    end
end


@testset "differentiation" begin
    @testset "gradient vs spline" begin

    end
end



