import DensityEstimators as DE
using QuadGK

@testset "Epanechnikov Kernel Tests" begin
    # Test case 1: u = 0.0 (center of the kernel)
    @test DE.kernel_epanechnikov(0.0) ≈ 0.75

    # Test case 2: u = 1.0 (boundary of the kernel)
    @test DE.kernel_epanechnikov(1.0) ≈ 0.0

    # Test case 3: u = -1.0 (boundary of the kernel)
    @test DE.kernel_epanechnikov(-1.0) ≈ 0.0

    # Test case 4: u = 0.5 (within the kernel)
    @test DE.kernel_epanechnikov(0.5) ≈ 0.75 * (1.0 - 0.5^2)

    # Test case 5: u = 1.5 (outside the kernel)
    @test DE.kernel_epanechnikov(1.5) ≈ 0.0

    # Test case 6: u = -1.5 (outside the kernel)
    @test DE.kernel_epanechnikov(-1.5) ≈ 0.0

    # Test case 7: u = 0.25 (within the kernel)
    @test DE.kernel_epanechnikov(0.25) ≈ 0.75 * (1.0 - 0.25^2)

    # Test case 8: u = -0.75 (within the kernel)
    @test DE.kernel_epanechnikov(-0.75) ≈ 0.75 * (1.0 - 0.75^2)
end



@testset "epanechnikov integral" begin
    k(a, b) = DE.integrate(DE.kernel_epanechnikov, a, b)
    k1(a, b) = quadgk(DE.kernel_epanechnikov, a, b)[1]


    # Test case 1: a = -1.0, b = 1.0
    @test k(-1.0, 1.0) ≈ k1(-1.0, 1.0) rtol=1e-2

    @test k(1.0, 2.5) ≈ 0 
    @test k(-2.5, -1) ≈ 0 

    for (a, b) in [(-0.9, -0.4), (-0.4, 0.2), (0.1, 0.34), (0.94, 0.99)]
        @test k(a, b) ≈ k1(a, b) rtol=1e-2
    end

    # test definition is normalized
    @test DE.normalize_kernel(DE.kernel_epanechnikov, 1.0) ≈ 1.0 rtol=1e-2
end


@testset "epanechnikov derivative" begin

end
