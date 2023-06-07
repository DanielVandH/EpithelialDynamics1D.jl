using EpithelialDynamics1D
using Test
using SafeTestsets

@safetestset "Problem" begin
    include("problem.jl")
end
@safetestset "Equations" begin
    include("equations.jl")
end
@safetestset "Proliferation" begin
    include("proliferation.jl")
end

@safetestset "Uniform" begin
    include("uniform.jl")
end
@safetestset "Step Function" begin
    include("step_function.jl")
end