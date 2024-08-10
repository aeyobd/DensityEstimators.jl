import StatsBase

@testset "bins_min_width_equal_number" begin
    # Test 1: Basic functionality
    @testset "Basic Functionality" begin
        x = [1.0, 2.0, 2.1, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        dx_min = 1.0
        N_per_bin_min = 2
        bins = Arya.bins_min_width_equal_number(x; dx_min=dx_min, N_per_bin_min=N_per_bin_min)
        @test length(bins) >= 1
        h = Arya.histogram(x, bins)
        @test all(h.values .>= N_per_bin_min)
        @test all(diff(bins) .>= dx_min)
    end

    # Test 2: Edge cases
    @testset "Edge Cases" begin
        # Empty array
        x = Float64[]
        dx_min = 1.0
        N_per_bin_min = 2
        bins = Arya.bins_min_width_equal_number(x; dx_min=dx_min, N_per_bin_min=N_per_bin_min)
        @test length(bins) == 0

        # Single element array
        x = [1.0]
        bins = Arya.bins_min_width_equal_number(x; dx_min=dx_min, N_per_bin_min=N_per_bin_min)
        @test length(bins) == 1

        # All elements the same
        x = [2.0, 2.0, 2.0]
        bins = Arya.bins_min_width_equal_number(x; dx_min=dx_min, N_per_bin_min=N_per_bin_min)
        @test length(bins) == 1
    end

    # Test 3: Parameter boundaries
    @testset "Parameter Boundaries" begin
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        dx_min = 0.0
        N_per_bin_min = 1
        bins = Arya.bins_min_width_equal_number(x; dx_min=dx_min, N_per_bin_min=N_per_bin_min)
        @test length(bins) == length(x) + 1

        dx_min = 5.0
        N_per_bin_min = 5
        bins = Arya.bins_min_width_equal_number(x; dx_min=dx_min, N_per_bin_min=N_per_bin_min)
        @test length(bins) == 2
        h = Arya.histogram(x, bins)
        @test all(h.values .>= 1)
    end

    # Test 4: Large array for performance
    @testset "Performance" begin
        x = randn(10000)
        dx_min = 0.1
        N_per_bin_min = 100
        @testset "Performance with large array" begin
            @time bins = Arya.bins_min_width_equal_number(x; dx_min=dx_min, N_per_bin_min=N_per_bin_min)
            @test length(bins) > 0
            h = Arya.histogram(x, bins)
            @test all(h.values .>= N_per_bin_min)
            @test all(diff(bins) .>= prevfloat(dx_min, 2))
        end
    end


    @testset "N equals zero" begin
        x = randn(10000)
        dx_min = 0.1
        N_per_bin_min = 0
        bins = Arya.bins_min_width_equal_number(x; dx_min=dx_min, N_per_bin_min=N_per_bin_min)
        @test length(bins) > 0
        h = Arya.histogram(x, bins)
        @test diff(bins[1:end-1]) ≈ fill(dx_min, length(bins)-2) rtol=1e-6
    end

    @testset "dx equals zero" begin
        x = randn(10000)
        dx_min = 0.0
        N_per_bin_min = 100
        bins = Arya.bins_min_width_equal_number(x; dx_min=dx_min, N_per_bin_min=N_per_bin_min)
        @test length(bins) > 0
        h = Arya.histogram(x, bins)
        @test all(h.values .>= N_per_bin_min)
    end
end
