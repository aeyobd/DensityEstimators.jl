import Random: randperm

@testset "knn" begin
    x = 0:0.1:1
    y = x .^ 2

    idxs, dists = Arya.knn(y, k=1)
    @test length(idxs) == length(x)
    @test length(dists) == length(x)

    N = length(x)
    idx_exp = vcat(2, 1:N-1)
    @test idxs == idx_exp
    dist_exp = vcat(0.01, diff(y))
    @test dists ≈ dist_exp

    perm = randperm(N)
    y = y[perm]
    idxs, dists = Arya.knn(y, k=1)
    @test dists ≈ dist_exp[perm]


    x = [1, 2, 2.5, -3.024, 7.23]

    idxs, dists = Arya.knn(x, k=2)
    @test length(idxs) == length(x)
    @test idxs[1] == 3

    idxs, dists = Arya.knn(x, k=3)
    @test idxs[1] == 4

    idxs, dists = Arya.knn(x, k=4)
    @test idxs[1] == 5
end
