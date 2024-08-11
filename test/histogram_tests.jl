import DensityEstimators: make_bins, bin_indices, bin_volumes, normalize, calc_hist_errors, simple_hist, histogram

@testset "bin_indicies" begin 
    @testset "bin index (left, default)" begin
        bins = [-Inf, -1.0, -0.5, 0.0, 0.00001, Inf]

        bin_index(x, bins) = DensityEstimators.bin_indices(x, bins)[1]

        @test bin_index(-1, bins) == 2
        @test bin_index(-1.1, bins) == 1
        @test bin_index(0.5, bins) == 5

        bins = [-1.2, -0.5, 0.0, 1e-5, 2.5]
        N = length(bins)

        @test bin_index(-2.0, bins) == 0
        @test bin_index(10., bins) == N 
        @test bin_index(NaN, bins) == N 
        @test bin_index(2.5, bins) == N - 1
        @test bin_index(-1.2, bins) == 1
        @test bin_index(-1., bins) == 1
        @test bin_index(1., bins) == 4
    end


    @testset "bin index right" begin
        bins = [-Inf, -1.0, -0.5, 0.0, 0.00001, Inf]

        bin_index(x, bins) = DensityEstimators.bin_indices(x, bins, :right)[1]
        @test bin_index(-1.0, bins) == 1
        @test bin_index(-1.1, bins) == 1
        @test bin_index(-0.9, bins) == 2
        @test bin_index(0.5, bins) == 5

        bins = [-1.2, -0.5, 0.0, 1e-5, 2.5]

        N = length(bins)
        @test bin_index(-2.0, bins) == 0
        @test bin_index(10.0, bins) == N
        @test bin_index(NaN, bins) == N
        @test bin_index(2.5, bins) == 4
        @test bin_index(-1.2, bins) == 1
        @test bin_index(-1.0, bins) == 1
        @test bin_index(1.0, bins) == 4
    end


    @testset "errors" begin
        @test_throws ArgumentError DensityEstimators.bin_indices([0], [1, 2, 3], :not_valid)
    end


    @testset "pathological" begin
        bins = []
        @test_throws ArgumentError DensityEstimators.bin_indices([0], bins)

        bins = [1]
        @test_throws ArgumentError DensityEstimators.bin_indices([0], bins)
    end
end


@testset "make_bins" begin
    @testset "number" begin
        x = []
        limits = (0, 3)
        bins = DensityEstimators.make_bins(x, limits, 3)
        @test bins ≈ [0, 1, 2, 3]

        limits = (-2, 1)
        bins = DensityEstimators.make_bins(x, limits, 6)
        @test bins ≈ -2:0.5:1 
        @test length(bins) == 7

    end


    @testset "bandwidth" begin
        x = []
        limits = (0, 1)
        bins = DensityEstimators.make_bins(x, limits, bandwidth=0.1)
        @test bins ≈ 0:0.1:1.1

        limits = (2, 2.313)
        bins = DensityEstimators.make_bins(x, limits, bandwidth=0.05)
        @test bins ≈ 2:0.05:2.363
    end


    @testset "function -> bandwidth" begin
        x = []
        limits = (0, 1)
        
        bw = 0.123
        f = x->0.123
        bins = DensityEstimators.make_bins(x, limits, bandwidth=f)
        @test bins ≈ 0:bw:(1+bw)


        limits = (-1.5, 1.5)
        bw = 0.2
        f = x->bw
        bins = DensityEstimators.make_bins(x, limits, bandwidth=f)
        @test bins ≈ -1.5:bw:(1.5+bw)
    end


    @testset "function -> bins" begin
        x = []
        limits = (0, 1)
        
        bw = 0.123
        bins = [0, 0.1, 0.2, 0.4, 1]
        f = x->bins
        bins = DensityEstimators.make_bins(x, limits, f)
        @test bins ≈ bins

    end


    @testset "array" begin
        x = [1,2]
        limits = (1, 2)
        
        bw = 0.123

        bins = make_bins(x, limits, [0, 0.1, 0.2, 0.4, 1])
    end


    @testset "exceptions" begin 
        x = [2, 5, 3, 4]
        limits = (0, 5)
        @test_throws ArgumentError make_bins(x, limits)
        @test_throws ArgumentError make_bins(x, limits; bandwidth="hehe")
        @test_throws ArgumentError make_bins(x, limits; bandwidth=x->"hehe")
        @test_throws ArgumentError make_bins(x, limits, x->"hehe")

        @test_throws ArgumentError make_bins(x, limits, [1,1,1,1,1])
        @test_throws ArgumentError make_bins(x, limits, [4,3,2,1])
        @test_throws ArgumentError make_bins(x, limits, [])
        @test_throws ArgumentError make_bins(x, limits, 0)
    end

    @testset "pathological" begin

    end
end


@testset "bin_volumes" begin
    @testset "simple cases" begin
        bins = [0, 1, 3, 5.5]
        volumes = DensityEstimators.bin_volumes(bins)
        @test volumes ≈ [1, 2, 2.5]
    end

    @testset "identities" begin
    end

    @testset "properties" begin
    end

    @testset "pathological" begin
        bins = Float64[]
        volumes = DensityEstimators.bin_volumes(bins)
        @test volumes == []


        bins = [0]
        volumes = DensityEstimators.bin_volumes(bins)
        @test volumes == []


        bins = [Inf, Inf, -Inf, NaN, NaN, 1., Inf]
        volumes = DensityEstimators.bin_volumes(bins)
        println(volumes)
        @test all(isequal.(volumes, [NaN, -Inf, NaN, NaN, NaN, Inf]))


        bins = [1, 1]
        volumes = DensityEstimators.bin_volumes(bins)
        @test volumes == [0]
    end
end


@testset "normalize" begin
    @testset "simple cases" begin
        hist = [1.0, 2, 0.2, 0]
        volumes = [1.0, 2.0, 0.5, π]

        hist_norm = DensityEstimators.normalize(hist, volumes, :none)
        @test hist_norm == hist

        hist_norm = DensityEstimators.normalize(hist, volumes, :pdf)
        @test hist_norm ≈ [1, 1, 0.4, 0] / 3.2

        hist_norm = DensityEstimators.normalize(hist, volumes, :density)
        @test hist_norm ≈ [1, 1, 0.4, 0]

        hist_norm = DensityEstimators.normalize(hist, volumes, :probabilitymass)
        @test hist_norm ≈ [1, 2, 0.2, 0] / 3.2
    end
        

    @testset "exceptions" begin
        hist = [1.0, 2.5, -0.5, 2.3]
        volumes = [1.0, 2.0, 0.5, 2]
        @test_throws ArgumentError DensityEstimators.normalize(hist, volumes, :gibberish)

        hist = [1.0, 2.5, -0.5, 2.3]
        volumes = [1.0, 2.0, 0.5]
        @test_throws DimensionMismatch DensityEstimators.normalize(hist, volumes, :pdf)
    end


    @testset "zero volume" begin
        # a little unphysical but worthwill to test
        
        hist = [1.0, 2.5, -0.5, 2.3]
        volumes = [1.0, 2.0, 0.5, 0]

        hist_norm = DensityEstimators.normalize(hist, volumes, :none)
        @test hist_norm == hist

        hist_norm = DensityEstimators.normalize(hist, volumes, :pdf)
        @test hist_norm ≈ [1, 1.25, -1, Inf] / 5.3

        hist_norm = DensityEstimators.normalize(hist, volumes, :density)
        @test hist_norm ≈ [1, 1.25, -1, Inf]

        hist_norm = DensityEstimators.normalize(hist, volumes, :probabilitymass)
        @test hist_norm ≈ [1, 2.5, -0.5, 2.3] / 5.3
    end

    @testset "infs" begin
        hist = [1.0, Inf, -0.5, 2.3]
        volumes = [1.0, 2.0, 0.5, 3.0]

        hist_norm = DensityEstimators.normalize(hist, volumes, :none)
        @test hist_norm == hist

        hist_norm = DensityEstimators.normalize(hist, volumes, :pdf)
        println(hist_norm)
        @test all(isequal(hist_norm, [0., NaN, -0., 0.]))

        hist_norm = DensityEstimators.normalize(hist, volumes, :density)
        @test hist_norm ≈ [1.0, Inf, -1.0, 2.3/3]

        hist_norm = DensityEstimators.normalize(hist, volumes, :probabilitymass)
        println(hist_norm)
        @test all(isequal.(hist_norm, [0., NaN, -0., 0.]))
    end


    @testset "nans" begin
        # a NaN should ruin everything
        hist = [1.0, NaN]
        volumes = [π, 2.0]

        hist_norm = DensityEstimators.normalize(hist, volumes, :none)
        @test all(hist_norm .=== hist)

        hist_norm = DensityEstimators.normalize(hist, volumes, :pdf)
        @test all(isnan.(hist_norm))

        hist_norm = DensityEstimators.normalize(hist, volumes, :density)
        @test hist_norm[1] ≈ 1/π
        @test hist_norm[2] === NaN

        hist_norm = DensityEstimators.normalize(hist, volumes, :probabilitymass)
        @test all(isnan.(hist_norm))
    end

    @testset "random: none" begin
        for _ in 1:100
            N = rand(1:100)
            hist = rand(N)
            volumes = rand(N)
            hist_norm = DensityEstimators.normalize(hist, volumes, :none)
            @test hist_norm == hist
        end
    end

    @testset "random: pdf" begin
        for _ in 1:100
            N = rand(1:100)
            hist = rand(N)
            volumes = rand(N)
            hist_norm = DensityEstimators.normalize(hist, volumes, :pdf)
            @test all(hist_norm .≥ 0)
            @test sum(hist_norm .* volumes) ≈ 1
        end
    end

    @testset "random: density" begin 
        for _ in 1:100 
            N = rand(1:100)
            hist = rand(N)
            volumes = rand(N)
            hist_norm = DensityEstimators.normalize(hist, volumes, :density)
            @test all(hist_norm .≥ 0)
            @test sum(hist_norm .* volumes) ≈ sum(hist)
        end
    end

    @testset "random: probabilitymass" begin
        for _ in 1:100
            N = rand(1:100)
            hist = rand(N)
            volumes = rand(N)
            hist_norm = DensityEstimators.normalize(hist, volumes, :probabilitymass)
            @test all(hist_norm .≥ 0)
            @test sum(hist_norm) ≈ 1
        end
    end


    @testset "pathological" begin
        hist = Real[]
        volumes = Real[]
        hist_norm = DensityEstimators.normalize(hist, volumes, :none)
        @test hist_norm == hist
    end
end



@testset "calc_hist_errors" begin
    @testset "poisson_errror" begin
        x = randn(1000)
        bins = -5:0.1:5
        hist, l, h = DensityEstimators.simple_hist(x, bins)

        hist_err = DensityEstimators.calc_hist_errors(x, bins, hist; errors = :poisson)

        @test all((hist_err .≥ 0) .| isnan.(hist_err))
        @test all((hist_err .≤ hist) .| isnan.(hist_err))

        filt = .! isnan.(hist_err)
        @test hist_err[filt] ≈ sqrt.(hist[filt])

        counts = copy(hist)

        hist = DensityEstimators.normalize(counts, DensityEstimators.bin_volumes(bins), :pdf)
        hist_err = DensityEstimators.calc_hist_errors(x, bins, hist; errors = :poisson)

        @test all((hist_err .≥ 0) .| isnan.(hist_err))
        @test all((hist_err .≤ hist) .| isnan.(hist_err))

        filt = .! isnan.(hist_err)
        @test hist_err[filt] ≈ hist[filt] ./ sqrt.(counts[filt])
    end

    @testset "no errors" begin
        x = randn(1000)
        bins = -5:0.1:5
        hist, l, h = DensityEstimators.simple_hist(x, bins)

        hist_err = DensityEstimators.calc_hist_errors(x, bins, hist; errors = :none)
        @test hist_err === nothing
    end

    @testset "exceptions" begin
        x = randn(10)
        bins = -1:0.1:1
        hist, _, _ = DensityEstimators.simple_hist(x, bins)

        @test_throws ArgumentError DensityEstimators.calc_hist_errors(x, bins, hist; errors = :gibberish)
    end
end



@testset "simple_hist" begin
    @testset "simple cases" begin
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        bins = [0, 5, 10]
        hist, l, h = DensityEstimators.simple_hist(x, bins; closed=:right)

        @test eltype(hist) <: Integer
        @test hist ≈ [5, 5]
    end

    @testset "properties" begin
        for i in 1:100
            N = rand(1:100)
            p = 10*randn()
            x = p * randn(N)
            bins = -10:0.1:10
            hist, l, h = DensityEstimators.simple_hist(x, bins)

            # mass conservation
            @test sum(hist) + l + h ≈ N
            @test l == sum(x .< bins[1])
            @test h == sum(x .≥ bins[end])

            
            # permutation invariance
            p = randperm(N)
            hist2, l2, h2 = DensityEstimators.simple_hist(x[p], bins)
            @test hist2 ≈ hist
            @test l2 == l
            @test h2 == h

            # nonnegativity
            @test all(hist .≥ 0)
        end
    end

    @testset "pathological" begin
        # TODO:
    end
end


@testset "simple_hist weighted" begin
    @testset "simple cases" begin
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        weights = [1, 2, 0, 0.5, 3, 1, 0, 0, 1, -1]

        bins = [0, 5, 10]
        hist, l, h = DensityEstimators.simple_hist(x, bins, weights; closed=:right)
        @test hist ≈ [6.5, 1]
        @test eltype(hist) <: AbstractFloat
    end

    @testset "identities" begin
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x = x[randperm(length(x))]
        weights = ones(length(x))
        bins = [0, 5, 10]
        hist, l, h = DensityEstimators.simple_hist(x, bins, weights; closed=:right)
        @test hist ≈ [5.0, 5.0]
    end

    @testset "properties" begin
        for i in 1:100
            N = rand(1:100)
            p = 10*randn()
            x = p * randn(N)
            w = rand(N)
            bins = -10:0.1:10

            hist, l, h = DensityEstimators.simple_hist(x, bins, w)

            # mass conservation
            @test sum(hist) + l + h ≈ sum(w)
            @test l ≈ sum(w[x .< bins[1]])
            @test h ≈ sum(w[x .≥ bins[end]])

            
            # permutation invariance
            p = randperm(N)
            hist2, l2, h2 = DensityEstimators.simple_hist(x[p], bins, w[p])
            @test hist2 ≈ hist
            @test l2 ≈ l
            @test h2 ≈ h

        end
    end

    @testset "exceptions" begin
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        weights = [1, 2, 0]

        bins = [0, 5, 10]
        @test_throws DimensionMismatch DensityEstimators.simple_hist(x, bins, weights; closed=:right)
    end

    @testset "pathological" begin
        # TODO:
    end
end


@testset "histogram" begin
    # these are more of integration tests
    #
    @testset "simple cases" begin
        x = [1, 2.1, 3.1, 3.2, 3.6, 4.92, 5.2, 8.2, 13.7]
        x = x[randperm(length(x))]
        bins = [2, 4, 5.5, 10]

        h = DensityEstimators.histogram(x, bins)

        @test h.bins ≈ bins
        @test h.values ≈ [4, 2, 1]
        @test h.err ≈ sqrt.([4, 2, 1])

        @test h.normalization == :none
        @test h.closed == :left

        h = DensityEstimators.histogram(x)
    end

    @testset "statistical" begin
        x = randn(10_000)
        bins = -3:0.2:3
        h = DensityEstimators.histogram(x, bins, normalization=:pdf)

        model(x) = exp(-x^2/2) / sqrt(2π)

        bin_mid = DensityEstimators.midpoints(h.bins)
        @test h.values ≈ model.(bin_mid) atol=0.1

        chi2 = sum((h.values .- model.(bin_mid)).^2 ./ h.err.^2)
        chi2 /= length(h.values)

        println(chi2)
        @test chi2  ≈ 1.0 atol=0.5
    end

end
