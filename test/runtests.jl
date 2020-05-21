using Thermodynamics
using Test

#--- TEST SET ------------------------------------------------------------------ #
function compute_pressure_ideal_gas_test()

    # create an ideal gass equation of state, compute pressure -
    R = 8.314e-2 # L bar mol^-1 K^-1
    eos_model = IdealGasEquationOfState(R)

    # set T and V -
    T = 293.15  # K
    V = 0.15    # L
    P_true = (R*T)/(V)
    P_calc = P(eos_model, V, T)

    # compute -
    return (P_calc == P_true)
end
#------------------------------------------------------------------------------- #

@testset "user_test_set" begin
    @test compute_pressure_ideal_gas_test() == true
end