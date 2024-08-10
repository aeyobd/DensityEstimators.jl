

@testset "bin index left" begin

    bins = [-Inf, -1.0, -0.5, 0.0, 0.00001, Inf]

    @test DensityEstimators.bin_index_left(bins, -1.0) == 2
    @test DensityEstimators.bin_index_left(bins, -1.1) == 1
    @test DensityEstimators.bin_index_left(bins, 0.5) == 5

    bins = [-1.2, -0.5, 0.0, 1e-5, 2.5]
    N = length(bins)

    @test DensityEstimators.bin_index_left(bins, -2.0) == 0
    @test DensityEstimators.bin_index_left(bins, 10) == N 
    @test DensityEstimators.bin_index_left(bins, NaN) == N 
    @test DensityEstimators.bin_index_left(bins, 2.5) == N - 1
    @test DensityEstimators.bin_index_left(bins, -1) == 1
    @test DensityEstimators.bin_index_left(bins, 1) == 4
end


@testset "bin index right" begin

    bins = [-Inf, -1.0, -0.5, 0.0, 0.00001, Inf]

    @test DensityEstimators.bin_index_right(bins, -1.0) == 1
    @test DensityEstimators.bin_index_right(bins, -1.1) == 1
    @test DensityEstimators.bin_index_right(bins, -0.9) == 2
    @test DensityEstimators.bin_index_right(bins, 0.5) == 5

    bins = [-1.2, -0.5, 0.0, 1e-5, 2.5]

    N = length(bins)
    @test DensityEstimators.bin_index_right(bins, -2.0) == 0
    @test DensityEstimators.bin_index_right(bins, 10) == N
    @test DensityEstimators.bin_index_right(bins, NaN) == N
    @test DensityEstimators.bin_index_right(bins, 2.5) == 4
    @test DensityEstimators.bin_index_right(bins, -1) == 1
    @test DensityEstimators.bin_index_right(bins, 1) == 4
end


@testset "make bins number" begin
    x = nothing
    limits = [0, 3]
    bins = DensityEstimators.make_bins(x, limits, 3)
    @test bins ≈ [0, 1, 2, 3]

    limits = [-2, 1]
    bins = DensityEstimators.make_bins(x, limits, 6)
    @test bins ≈ -2:0.5:1 
    @test length(bins) == 7

end


@testset "make bins bandwidth" begin
    x = nothing
    limits = [0, 1]
    bins = DensityEstimators.make_bins(x, limits, bandwidth=0.1)
    @test bins ≈ 0:0.1:1.1

    limits = [2, 2.313]
    bins = DensityEstimators.make_bins(x, limits, bandwidth=0.05)
    @test bins ≈ 2:0.05:2.363
end


@testset "make bins function -> bandwidth" begin
    x = nothing
    limits = [0, 1]
    
    bw = 0.123
    f = x->0.123
    bins = DensityEstimators.make_bins(x, limits, f)
    @test bins ≈ 0:bw:(1+bw)


    limits = [-1.5, 1.5]
    bw = 0.2
    f = x->bw
    bins = DensityEstimators.make_bins(x, limits, f)
    @test bins ≈ -1.5:bw:(1.5+bw)
end


@testset "make bins function -> bins" begin
    x = nothing
    limits = [0, 1]
    
    bw = 0.123
    bins = [0, 0.1, 0.2, 0.4, 1]
    f = x->bins
    bins = DensityEstimators.make_bins(x, limits, f)
    @test bins ≈ bins

end


@testset "make bins array" begin
    x = nothing
    limits = [0, 1]
    
    bw = 0.123
end

@testset "midpoints" begin
    x = 1:5
    @test DensityEstimators.midpoints(x) ≈ 1.5:1:4.5

    x = [1, -2, 4]
    @test DensityEstimators.midpoints(x) ≈ [-0.5, 1]
end
