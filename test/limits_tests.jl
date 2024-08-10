
x = [-0.2, 0.1, -0.546, π, 3, 1.56, 20.75, NaN, Inf, -Inf]
x_min = -0.546
x_max = 20.75

@testset "no limits" begin
    limits = Arya.calc_limits(x, filt_inf=true)
    @test limits[1] == x_min
    @test limits[2] == x_max


    limits = Arya.calc_limits(x, nothing, filt_inf=true)
    @test limits[1] == x_min
    @test limits[2] == x_max


    limits = Arya.calc_limits(x, (nothing, nothing), filt_inf=true)
    @test limits[1] == x_min
    @test limits[2] == x_max
end

@testset "only upper" begin
    for up in [1.23, 2, 50]
        limits = Arya.calc_limits(x, (nothing, up), filt_inf=true)
        @test limits[1] == x_min
        @test limits[2] === up
        println(limits)
    end

    for up in [1.23, 2, 50, Inf]
        limits = Arya.calc_limits(x, (nothing, up), filt_inf=false)
        @test limits[1] == -Inf
        @test limits[2] === min(up, Inf)
    end
end


@testset "only lower" begin
    for lower in [-1, 2, -20.23]
        limits = Arya.calc_limits(x, (lower, nothing), filt_inf=true)
        @test limits[1] === lower
        @test limits[2] == x_max
    end
end




y = [0.2, 0.4, -0.2, 2π, -4, 1.56, 9.75, NaN, Inf, -Inf]
y_min = -4
y_max = 9.75

@testset "both limits" begin
    xlims = Arya.calc_limits(x, filt_inf=true)
    ylims = Arya.calc_limits(y, filt_inf=true)
    @test xlims[1] == x_min
    @test xlims[2] == x_max
    @test ylims[1] == y_min
    @test ylims[2] == y_max

end
