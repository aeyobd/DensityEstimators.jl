
include("setup.jl")


tests = ["limits", "bandwidth", "histogram", "interface", "bayesian_blocks", "histogram2d",
         "kde2d",
        ]

for test in tests
    @testset "$test" begin
        include("$(test)_tests.jl")
    end
end
