function obj_function_simple_vstar_lk(x, eos::LeeKeslerEquationOfState, P::Float64, T::Float64)

    # get x values -
    Vstar = x[1]

    # get parameters from the model -
    parameter_dict = eos.parameters
    R = parameter_dict["R"]
    Tc = parameter_dict["Tc"]
    Pc = parameter_dict["Pc"]
    Vc = parameter_dict["Vc"]
    ω = parameter_dict["ω"]

    # compute reduced T -
    Tr = T/Tc
    Pr = P/Pc

    # compute simple fluid model Z -
    b1 = parameter_dict["simple_fluid_parameters"]["b1"]
    b2 = parameter_dict["simple_fluid_parameters"]["b2"]
    b3 = parameter_dict["simple_fluid_parameters"]["b3"]
    b4 = parameter_dict["simple_fluid_parameters"]["b4"]
    c1 = parameter_dict["simple_fluid_parameters"]["c1"]
    c2 = parameter_dict["simple_fluid_parameters"]["c2"]
    c3 = parameter_dict["simple_fluid_parameters"]["c3"]
    c4 = parameter_dict["simple_fluid_parameters"]["c4"]
    d1 = parameter_dict["simple_fluid_parameters"]["d1"]
    d2 = parameter_dict["simple_fluid_parameters"]["d2"]
    β = parameter_dict["simple_fluid_parameters"]["β"]
    𝛄 = parameter_dict["simple_fluid_parameters"]["𝛄"]

    B = b1 - b2/Tr - b3/(Tr^2) - b4/(Tr^3)
    C = c1 - c2/Tr + c3/(Tr^3)
    D = d1 + d2/Tr
    err = (Vstar*Pr)/(Tr) - (1 + B/Vstar+C/(Vstar^2)+D/(Vstar^5) + (c4/((Tr^3)*(Vstar^2)))*(β+𝛄/(Vstar^2))*exp(-𝛄/(Vstar^2))) 
    
    # penalty terms -
    p_term_1 = max(0,-1*Vstar)

    # return -
    return (err*err) + 100000.0*(p_term_1)
end

function obj_function_complex_vstar_lk(x, eos::LeeKeslerEquationOfState, P::Float64, T::Float64)

    # get x values -
    Vstar = x[1]

    # get parameters from the model -
    parameter_dict = eos.parameters
    R = parameter_dict["R"]
    Tc = parameter_dict["Tc"]
    Pc = parameter_dict["Pc"]
    Vc = parameter_dict["Vc"]
    ω = parameter_dict["ω"]

    # compute reduced T -
    Tr = T/Tc
    Pr = P/Pc

    # compute simple fluid model Z -
    b1 = parameter_dict["reference_fluid_parameters"]["b1"]
    b2 = parameter_dict["reference_fluid_parameters"]["b2"]
    b3 = parameter_dict["reference_fluid_parameters"]["b3"]
    b4 = parameter_dict["reference_fluid_parameters"]["b4"]
    c1 = parameter_dict["reference_fluid_parameters"]["c1"]
    c2 = parameter_dict["reference_fluid_parameters"]["c2"]
    c3 = parameter_dict["reference_fluid_parameters"]["c3"]
    c4 = parameter_dict["reference_fluid_parameters"]["c4"]
    d1 = parameter_dict["reference_fluid_parameters"]["d1"]
    d2 = parameter_dict["reference_fluid_parameters"]["d2"]
    β = parameter_dict["reference_fluid_parameters"]["β"]
    𝛄 = parameter_dict["reference_fluid_parameters"]["𝛄"]

    B = b1 - b2/Tr - b3/(Tr^2) - b4/(Tr^3)
    C = c1 - c2/Tr + c3/(Tr^3)
    D = d1 + d2/Tr
    err = (Vstar*Pr)/(Tr) - (1 + B/Vstar+C/(Vstar^2)+D/(Vstar^5) + (c4/((Tr^3)*(Vstar^2)))*(β+𝛄/(Vstar^2))*exp(-𝛄/(Vstar^2))) 
    
    # penalty terms -
    p_term_1 = max(0,-1*Vstar)

    # return -
    return (err*err) + 100000.0*(p_term_1)
end

function compute_simple_compressibility_lk(eos::LeeKeslerEquationOfState, P::Float64, T::Float64; V::Float64=0.01)

    # get parameters -
    parameter_dict = eos.parameters
    R = parameter_dict["R"]
    Tc = parameter_dict["Tc"]
    Pc = parameter_dict["Pc"]
    Vc = parameter_dict["Vc"]
    ω = parameter_dict["ω"]
 
    # setup the calculation -
    xinitial = [V]
    OF(p) = obj_function_simple_vstar_lk(p, eos, P, T)
     
    # call the optimizer -
    opt_result = optimize(OF,xinitial,BFGS())
 
    # what is Vstar?
    Vstar = Optim.minimizer(opt_result)[1]
    
    # compute Z -
    Tr = T/Tc
    Pr = P/Pc
    Z0 = (Vstar*Pr)/(Tr)

    # return -
    return Z0
end

function compute_complex_compressibility_lk(eos::LeeKeslerEquationOfState, Z0::Float64,P::Float64,T::Float64; V::Float64 = 0.01)

    # get parameters -
    parameter_dict = eos.parameters
    R = parameter_dict["R"]
    Tc = parameter_dict["Tc"]
    Pc = parameter_dict["Pc"]
    Vc = parameter_dict["Vc"]
    ω = parameter_dict["ω"]
    ωr = 0.3978
 
    # setup the calculation -
    xinitial = [V]
    OF(p) = obj_function_complex_vstar_lk(p, eos, P, T)
     
    # call the optimizer -
    opt_result = optimize(OF,xinitial,BFGS())
 
    # what is Vstar?
    Vstar = Optim.minimizer(opt_result)[1]

    # compute Z -
    Tr = T/Tc
    Pr = P/Pc
    Zr = (Vstar*Pr)/(Tr)

    # compute Z1 -
    Z1 = (1/ωr)*(Zr - Z0)

    # return -
    return Z1
end