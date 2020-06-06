function buildVanDerWaalsEquationOfState(fluid::SingleComponentWorkingFluid)::VanDerWaalsEquationOfState

    # compute the analytical vdw parameters -
    # get the critical parameters -
    Tc = fluid.Tc
    Pc = fluid.Pc
    Vc = fluid.Vc
    R = fluid.R

    # compute a, b from the critical parameters -
    a = (27/64)*((R^2)*(Tc^2))/(Pc)
    b = (1/8)*(R)*(Tc/Pc)

    # build a vdw EOS model -
    eos_model = VanDerWaalsEquationOfState(R,a,b)

    # return -
    return eos_model
end

function buildLeeKeslerEquationOfState(fluid::SingleComponentWorkingFluid)::LeeKeslerEquationOfState

    # get the critcal properties, and other stuff -
    Tc = fluid.Tc
    Pc = fluid.Pc
    Vc = fluid.Vc
    ω = fluid.ω
    R = fluid.R

end

function buildPengRobinsonEquationOfState(fluid::SingleComponentWorkingFluid)::PengRobinsonEquationOfState

    # get the critcal properties, and other stuff -
    Tc = fluid.Tc
    Pc = fluid.Pc
    Vc = fluid.Vc
    ω = fluid.ω
    R = fluid.R

    # compute a, b and kappa -
    a = 0.45724*((R^2)*(Tc^2))/(Pc)
    b = 0.07780*((R)*Tc)/(Pc)
    𝝹 = 0.37464 + 1.54226*(ω) - 0.26992*(ω)^2
    

end