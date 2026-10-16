using TestItemRunner

@testmodule Reference begin
    include(joinpath(@__DIR__, "reference.jl"))
end

@run_package_tests
