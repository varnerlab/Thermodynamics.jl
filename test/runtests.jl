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
    P_calc = pressure(eos_model, V, T)

    # compute -
    return (P_calc == P_true)
end

function compute_pressure_vdw_test()

    # create an ideal gass equation of state, compute pressure -
    R = 8.314e-2 # L bar mol^-1 K^-1
    a = 6.343    # bar L^2 mol^-2
    b = 0.05422  # L mol^-1
    eos_model = VanDerWaalsEquationOfState(R,a,b)

    # compute real pressure -
    T = 293.15  # K
    V = 0.15    # L
    P_true = (R*T)/(V-b) - a/V^2
    P_calc = pressure(eos_model,V,T)

    # compute -
    return (P_calc == P_true)
end
#------------------------------------------------------------------------------- #

@testset "user_test_set" begin
    @test compute_pressure_ideal_gas_test() == true
    @test compute_pressure_vdw_test() == true
end