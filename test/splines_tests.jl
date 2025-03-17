import DensityEstimators as DE
using QuadGK: quadgk

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

    @testset "order 0" begin
        t = [-56.2, -0.8, 20.4, 100.0]
        α = [1.0, 0.0, 3.0, -1.5]
        k = 0
        @test_throws ArgumentError DE.BSpline(t, α, k)
    end

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


@testset "differentiation" begin
    @testset "gradient vs spline" begin
        s = DE.BSpline([0.1, 0.35, 2.5, 7.2, 9.0, 12.3], [-0.2, 0.9, 0.4, 0.11], 2)
        ds = DE.derivative(s)

        x = LinRange(0, 12.57, 100)

        dsx = ds.(x)
        h = 1e-8
        dsx1 = [(s(x[i] + h) - s(x[i])) / h for i in eachindex(x)]

        @test dsx ≈ dsx1 rtol=1e-6


        s = DE.BSpline([0.1, 0.35, 2.5, 3.4, 7.2, 9.0, 12.3], [-0.2, 0.9, 0.4, 0.11], 3)
        ds = DE.derivative(s)

        x = LinRange(0, 12.57, 100)

        dsx = ds.(x)
        h = 1e-8
        dsx1 = [(s(x[i] + h) - s(x[i])) / h for i in eachindex(x)]

        @test dsx ≈ dsx1 rtol=1e-6
    end
end


@testset "integration" begin

    @testset "quadgk" begin
        s = DE.BSpline([0.0, 0.2, 0.5, 1.2, 1.8, 2.5, 5.1], [2.2, 1.05, -0.28, 1.3], 3)
        si = DE.integral(s)

        x1 = LinRange(0.0, 4.9, 10)
        x2 = x1 .+ LinRange(0.2, 0.1, 10)

        @test [quadgk(s, 0, x)[1] for x in x1] ≈ [si(x) for x in x1] rtol=1e-8

        @test [quadgk(s, x1[i], x2[i])[1] for i in eachindex(x1)] ≈ [si(x2[i]) - si(x1[i]) for i in eachindex(x1)] rtol=1e-8


        @test DE.area_of(s) ≈ si(5.1) rtol=1e-8
    end
end





