using Test

import Arya: histogram2d, Histogram2D

# Assuming the Histogram2D and histogram2d definitions are already imported.

# Helper function to compare two Histogram2D structs
function compare_histograms(h1::Histogram2D, h2::Histogram2D)
    @test h1.xbins == h2.xbins
    @test h1.ybins == h2.ybins
    @test h1.values == h2.values
end

# Test 1: Simple 2x2 histogram without weights and limits
@testset "test_simple_histogram" begin
    x = [1.0, 2.0, 1.0, 2.0]
    y = [1.0, 1.0, 2.0, 2.0]
    bins = ([0.5, 1.5, 2.5], [0.5, 1.5, 2.5])
    h_expected = Histogram2D(
        xbins = [0.5, 1.5, 2.5],
        ybins = [0.5, 1.5, 2.5],
        values = [1 1; 1 1]
    )
    h = histogram2d(x, y, bins)
    compare_histograms(h, h_expected)
end

# Test 2: 2x2 histogram with weights
@testset "test_histogram_with_weights" begin
    x = [1.0, 2.0, 1.0, 2.0]
    y = [1.0, 1.0, 2.0, 2.0]
    bins = ([0.5, 1.5, 2.5], [0.5, 1.5, 2.5])
    weights = [1.0, 0.5, 2.0, 1.5]
    h_expected = Histogram2D(
        xbins = [0.5, 
                 1.5, 
                 2.5],
        ybins = [0.5, 1.5, 2.5],
        values = [1.0   2.0
                  0.5   1.5]
    )
    h = histogram2d(x, y, bins, weights=weights)
    compare_histograms(h, h_expected)
end

# Test 3: Histogram with limits
@testset "test_histogram_with_limits" begin
    x = [1.0, 2.0, 1.0, 2.0, 3.0, 4.0]
    y = [1.0, 1.0, 2.0, 2.0, 3.0, 4.0]
    bins = 2
    limits = ((1.0, 2.0), (1.0, 2.0))
    h_expected = Histogram2D(
        xbins = [1.0, 1.5, 2.0],
        ybins = [1.0, 1.5, 2.0],
        values = [1 1; 1 1]
    )
    h = histogram2d(x, y, bins; limits=limits)
    compare_histograms(h, h_expected)
end

# Test 4: Histogram with bin edges as a vector
@testset "test_histogram_with_vector_bins" begin
    x = [1.0, 2.0, 1.0, 2.0]
    y = [1.0, 1.0, 2.0, 2.0]
    bins = [0.5, 1.5, 2.5]
    h_expected = Histogram2D(
        xbins = [0.5, 1.5, 2.5],
        ybins = [0.5, 1.5, 2.5],
        values = [1 1; 1 1]
    )
    h = histogram2d(x, y, bins)
    compare_histograms(h, h_expected)
end

# Test 5: Histogram with infinite values
@testset "test_histogram_with_infinite_values" begin
    x = [1.0, 
         2.0, 
         Inf, 
         -Inf]
    y = [1.0, 1.0, Inf, -Inf]
    bins = ([0.5, 1.5, 2.5], [0.5, 1.5, 2.5])
    h_expected = Histogram2D(
        xbins = [0.5, 1.5, 2.5],
        ybins = [0.5, 1.5, 2.5],
        values = [1 0
                  1 0]
    )
    h = histogram2d(x, y, bins)
    compare_histograms(h, h_expected)
end



@testset "test_histogram_normalization" begin

end
