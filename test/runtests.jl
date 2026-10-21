using Optimizer
using Test

@testset "Optimizer.jl" begin

    @time display(brentqv(sin, -4π:0.1:4π, 10) )
    @time display(brentqv(x -> x^2 - 1, -2:0.01:2, 10) )
end
