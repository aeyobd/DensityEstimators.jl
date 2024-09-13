@testset "binned statistics simple" begin

    bins = [0,  2,                    3.4,                          6,10,12]
    x = [1.0, 2.0, 2.7, 3.0, 3.1, 3.23, 3.45, 3.67, 5.10, 10]
    y = [0.0, 1.0, 0.5, -1., 2.0, 7.00, 0.20, -2.0, -0.1, π ]
    w = [0.0, 0.0, 0.2,  .1, 1.0, 0.00, 0.20, 0.05, -0.1, √2 ]

    expected = [0.5, 2.125, -1.9/3, π, NaN]

    idx = randperm(length(x))
    binned_stats = DensityEstimators.binned_statistic(x[idx], y[idx], bins, statistic=DensityEstimators.mean, closed=:right)
    @test binned_stats[1:end-1] ≈ expected[1:end-1]
    @test isnan(binned_stats[end])

    binned_stats = DensityEstimators.binned_statistic(x[idx], y[idx], bins, statistic=DensityEstimators.median, closed=:right)

    expected = [0.5, 1.25, -0.1, π, NaN]
    @test binned_stats[1:end-1] ≈ expected[1:end-1]
    @test isnan(binned_stats[end])

    binned_stats = DensityEstimators.binned_statistic(x[idx], y[idx], bins, statistic=DensityEstimators.mean, closed=:right, weights=w[idx])

    expected = [NaN, 2/1.3, -0.05/0.15, π, NaN]
    @test isnan(binned_stats[1])
    @test isnan(binned_stats[end])
    @test binned_stats[2:end-1] ≈ expected[2:end-1]
end


@testset "binned_statistic pathological" begin
    #@test false broken=true
end

@testset "binned_statistic errors" begin
    #@test false broken=true
end


@testset "binned statistic 2d" begin
    @testset "simple" begin
        x = [1.0, 2.0, 2.7, 3.0, 3.1, 3.23, 3.45, 3.67, 5.10, 10]
        y = [-1.0, 2, 2.5, -0.1, 0.5, 0.7,  0.9,  1.0,  1.1,  1.2]
        # bin indicies
        #    [1,   1,   2,   2,   2,   2,   3,   3,   3,  4]
        #    [1,   3,   4,   1,   2,   2,   2,   2,   3,  3]
        z = [0.0, 1.0, 0.5, -1., 2.0, 7.00, 0.20, -2.0, -0.1, π ]
        w = [0.0, 0.2, 0.2,  .1, 1.0, 0.00, 0.20, 0.05, -0.1, √2 ]

        bins = ([0, 2, 3.4, 6, 10], [-2, 0, 1, 2, 3])

        expected = [0   NaN  1.0  NaN
                    -1  4.5  NaN  0.5
                    NaN -0.9 -0.1 NaN
                    NaN NaN  π    NaN]



        idx = randperm(length(x))
        x = x[idx]
        y = y[idx]
        z = z[idx]
        w = w[idx]

        binned_stats = DensityEstimators.binned_statistic_2d(x, y, z, bins, statistic=DensityEstimators.mean, closed=:right)

        println("expected = $expected")
        println("binned_stats = $binned_stats")

        nan_idx = isnan.(expected)
        @test binned_stats[.!nan_idx] ≈ expected[.!nan_idx]
        @test all(isnan.(binned_stats[nan_idx]))

        binned_stats = DensityEstimators.binned_statistic_2d(x, y, z, bins, statistic=DensityEstimators.mean, closed=:right, weights=w)

        expected = [NaN NaN   1.0  NaN
                    -1  2.0   NaN  0.5
                    NaN -0.24 -0.1 NaN
                    NaN NaN    π   NaN]

        nan_idx = isnan.(expected)
        @test binned_stats[.!nan_idx] ≈ expected[.!nan_idx]
        @test all(isnan.(binned_stats[nan_idx]))

    end

end
