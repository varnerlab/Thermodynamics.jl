function buildSingleComponentWorkingFluid(propertyDictionary::Dict{AbstractString,Any})::SingleComponentWorkingFluid

    # get stuff from dictionary -
    Tc = propertyDictionary["critical_temperature"]
    Pc = propertyDictionary["critical_pressure"]
    Vc = propertyDictionary["critical_volume"]
    ω = propertyDictionary["accentric_factor"]
    R = propertyDictionary["ideal_gas_constant"]

    # build -
    return SingleComponentWorkingFluid(Tc,Pc,Vc,ω;R=R)
end


function buildRedlichKwongEquationOfState(fluid::SingleComponentWorkingFluid)::RedlichKwongEquationOfState

    # compute the analytical RK parameters -
    # get the critical parameters -
    Tc = fluid.Tc
    Pc = fluid.Pc
    Vc = fluid.Vc
    R = fluid.R

    # compute a, b from the critical parameters -
    a = 0.42748*((R^2)*(Tc)^2.5)/(Pc)
    b = 0.08664*(R*Tc)/Pc

    # build a RK EOS model -
    eos_model = RedlichKwongEquationOfState(R,a,b)

    # return -
    return eos_model
end

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

    # setup parameters -
    # Table 1 LK AIChE 21:511 1975
    # simple parameters -
    simple_fluid_parameters = Dict{String,Any}()
    simple_fluid_parameters["b1"] = 0.1181193
    simple_fluid_parameters["b2"] = 0.265728
    simple_fluid_parameters["b3"] = 0.154790
    simple_fluid_parameters["b4"] = 0.030323
    simple_fluid_parameters["c1"] = 0.0236744
    simple_fluid_parameters["c2"] = 0.0186984
    simple_fluid_parameters["c3"] = 0.0
    simple_fluid_parameters["c4"] = 0.042724
    simple_fluid_parameters["d1"] = 0.155488e-4
    simple_fluid_parameters["d2"] = 0.623689e-4
    simple_fluid_parameters["β"] = 0.65392
    simple_fluid_parameters["𝛄"] = 0.060167

    # reference parameters -
    reference_fluid_parameters = Dict{String,Any}()
    reference_fluid_parameters["b1"] = 0.2026579
    reference_fluid_parameters["b2"] = 0.331511
    reference_fluid_parameters["b3"] = 0.027655
    reference_fluid_parameters["b4"] = 0.203488
    reference_fluid_parameters["c1"] = 0.0313385
    reference_fluid_parameters["c2"] = 0.0503618
    reference_fluid_parameters["c3"] = 0.016901
    reference_fluid_parameters["c4"] = 0.041577
    reference_fluid_parameters["d1"] = 0.48736e-4
    reference_fluid_parameters["d2"] = 0.0740336e-4
    reference_fluid_parameters["β"] = 1.226
    reference_fluid_parameters["𝛄"] = 0.03754
    reference_fluid_parameters["Tc"] = 568.8
    reference_fluid_parameters["Pc"] = 2.49
    reference_fluid_parameters["Vc"] = 0.492
    reference_fluid_parameters["ωr"] = 0.3978

    # package the parameters -
    model_parameters = Dict{String,Any}()
    model_parameters["simple_fluid_parameters"] = simple_fluid_parameters
    model_parameters["reference_fluid_parameters"] = reference_fluid_parameters
    model_parameters["R"] = R
    model_parameters["Tc"] = Tc
    model_parameters["Pc"] = Pc
    model_parameters["Vc"] = Vc
    model_parameters["ω"] = ω

    # return -
    return LeeKeslerEquationOfState(model_parameters)
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

    # build a vdw EOS model -
    eos_model = PengRobinsonEquationOfState(R,a,b,𝝹,Tc)

    # return -
    return eos_model
end

function buildMartinHouEquationOfState(fluid::SingleComponentWorkingFluid)::MartinHouEquationOfState
end