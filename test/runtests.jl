include("setup.jl")
using Random

Random.seed!(314159)


tests = ["limits", "bandwidth", "histogram", "interface", "bayesian_blocks", "histogram2d",
         "kde2d", "binned_statistics", "splines"
        ]

for test in tests
    @testset "$test" begin
        include("$(test)_tests.jl")
    end
end
