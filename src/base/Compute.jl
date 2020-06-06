# --- PRIVATE METHODS ------------------------------------------------------ #
function obj_function_volume_vdw(x, eos::VanDerWaalsEquationOfState, P::Float64, T::Float64)

    # get R, a and b from the eos wrapper -
    R = eos.R
    a = eos.a
    b = eos.b

    # get x values -
    V = x[1]

    # compute terms -
    repulsion = (R*T)/(V-b)
    attraction = a/(V^2)

    # compute error -
    err = (P - repulsion + attraction)

    # penalty terms -
    p_term_1 = max(0,-1*V)
    p_term_2 = max(0,-1*(V-b))

    # return -
    return (err*err) + 10000.0*(p_term_1 + p_term_2)
end
# -------------------------------------------------------------------------- #

# --- Ideal gas EOS ------------------------------------------------------------------------------- #
function compute_volume_ideal_gas(eos::IdealGasEquationOfState, P::Float64, T::Float64)::Float64

    # check - is V = 0
    if (V==0.0)
        throw(ArgumentError("Volume cannot be zero!"))
    end

    # get R from the eos wrapper -
    R = eos.R

    # compute the pressure -
    return (R/P)*T
end

function compute_pressure_ideal_gas(eos::IdealGasEquationOfState, V::Float64, T::Float64)::Float64

    # check - is V = 0
    if (V==0.0)
        throw(ArgumentError("Volume cannot be zero!"))
    end

    # get R from the eos wrapper -
    R = eos.R

    # compute the pressure -
    return (R/V)*T
end

function compute_temperature_ideal_gas(eos::IdealGasEquationOfState, P::Float64, V::Float64)::Float64

    # check - is V = 0
    if (V==0.0)
        throw(ArgumentError("Volume cannot be zero!"))
    end

    # get R from the eos wrapper -
    R = eos.R

    # compute the pressure -
    return (P*V)/(R)
end

# --- Van der Waals EOS ----------------------------------------------------------------------------- #
function compute_volume_vdw(eos::VanDerWaalsEquationOfState, P::Float64, T::Float64; V::Float64 = 1.0)::Float64

    # setup the calculation -
    xinitial = [V]
    OF(p) = obj_function_volume_vdw(p, eos, P, T)
    
    # call the optimizer -
    opt_result = optimize(OF,xinitial,BFGS())
    
    # return -
    return Optim.minimizer(opt_result)[1]
end

function compute_pressure_vdw(eos::VanDerWaalsEquationOfState, V::Float64, T::Float64)::Float64

    # get R, a and b from the eos wrapper -
    R = eos.R
    a = eos.a
    b = eos.b

    # check - is V = 0
    if (V==0.0)
        throw(ArgumentError("Volume == 0: attraction term is singular"))
    elseif ((V-b) == 0.0)
        throw(ArgumentError("Volume == b: singular repulsion term"))
    end

    # compute terms -
    repulsion = (R*T)/(V-b)
    attraction = a/(V^2)

    # return pressure -
    return (repulsion - attraction)
end

function compute_temperature_vdw(eos::VanDerWaalsEquationOfState, P::Float64, V::Float64)::Float64

    # get R, a and b from the eos wrapper -
    R = eos.R
    a = eos.a
    b = eos.b

    # check - is V = 0
    if (V==0.0)
        throw(ArgumentError("Volume == 0: attraction term is singular"))
    elseif ((V-b) == 0.0)
        throw(ArgumentError("Volume == b: singular repulsion term"))
    end

    # return temperature -
    return ((V-b)/R)*(P+a/V^2)
end

# --- Peng Robinson EOS --------------------------------------------------------------------------- #
function compute_pressure_pr(eos::PengRobinsonEquationOfState, V::Float64, T::Float64)::Float64

    # get parameters from eos -
    R = eos.R
    a = eos.a
    b = eos.b
    𝞳 = eos.𝞳
    Tc = eos.Tc

    # check - is V = 0, or V - b = 0?
    if (V==0.0)
        throw(ArgumentError("Volume == 0: attraction term is singular"))
    elseif ((V-b) == 0.0)
        throw(ArgumentError("Volume == b: singular repulsion term"))
    end

    # compute terms -
    Tr = T/Tc
    𝜶 = (1+𝞳*(1-Tr^0.5))^2
    repulsion = (R*T)/(V-b)
    attraction = (a*𝜶)/(V^2+2*b*V-b^2)

    # return -
    return (repulsion - attraction)
end

function compute_volume_pr(eos::PengRobinsonEquationOfState, P::Float64, T::Float64)::Float64
end

function compute_temperature_pr(eos::PengRobinsonEquationOfState, P::Float64, V::Float64)::Float64
end

# --- Martin-Hou EOS ----------------------------------------------------------------------------- #
function compute_volume_mh(eos::MartinHouEquationOfState, V::Float64, T::Float64)::Float64
    
    # get parameters from the model -
    parameter_dict = eos.parameters
    R = parameter_dict["R"]
    A2 = parameter_dict["A2"]
    B2 = parameter_dict["B2"]
    C2 = parameter_dict["C2"]
    A3 = parameter_dict["A3"]
    B3 = parameter_dict["B3"]
    C3 = parameter_dict["C3"]
    A4 = parameter_dict["A4"]
    B4 = parameter_dict["B4"]
    C4 = parameter_dict["C4"]
    A5 = parameter_dict["A5"]
    B5 = parameter_dict["B5"]
    C5 = parameter_dict["C5"]
    b = parameter_dict["b"]
    k = parameter_dict["k"]
    Tc = parameter_dict["Tc"]

    # compute reduced T -
    Tr = T/Tc

    # check - is V = 0
    if (V==0.0)
        throw(ArgumentError("Volume == 0: attraction term is singular"))
    elseif ((V-b) == 0.0)
        throw(ArgumentError("Volume == b: singular repulsion term"))
    end

    # compute terms -
    term_array = Float64[]
    term_1 = (R*T)/(V-b)
    term_2 = (A2+B2*T+C2*exp(-k*Tr))/((V-b)^2)
    term_3 = (A3+B3*T+C3*exp(-k*Tr))/((V-b)^3)
    term_4 = (A4+B4*T+C4*exp(-k*Tr))/((V-b)^4)
    term_5 = (A5+B5*T+C5*exp(-k*Tr))/((V-b)^5)
    
    # cache -
    push!(term_array, term_1)
    push!(term_array, term_2)
    push!(term_array, term_3)
    push!(term_array, term_4)
    push!(term_array, term_5)

    # return -
    return sum(term_array)
end

function compute_pressure_mh(eos::MartinHouEquationOfState, V::Float64, T::Float64)::Float64
    
    # get parameters from the model -
    parameter_dict = eos.parameters
    R = parameter_dict["R"]
    A2 = parameter_dict["A2"]
    B2 = parameter_dict["B2"]
    C2 = parameter_dict["C2"]
    A3 = parameter_dict["A3"]
    B3 = parameter_dict["B3"]
    C3 = parameter_dict["C3"]
    A4 = parameter_dict["A4"]
    B4 = parameter_dict["B4"]
    C4 = parameter_dict["C4"]
    A5 = parameter_dict["A5"]
    B5 = parameter_dict["B5"]
    C5 = parameter_dict["C5"]
    b = parameter_dict["b"]
    k = parameter_dict["k"]
    Tc = parameter_dict["Tc"]

    # compute reduced T -
    Tr = T/Tc

    # check - is V = 0
    if (V==0.0)
        throw(ArgumentError("Volume == 0: attraction term is singular"))
    elseif ((V-b) == 0.0)
        throw(ArgumentError("Volume == b: singular repulsion term"))
    end

    # compute terms -
    term_array = Float64[]
    term_1 = (R*T)/(V-b)
    term_2 = (A2+B2*T+C2*exp(-k*Tr))/((V-b)^2)
    term_3 = (A3+B3*T+C3*exp(-k*Tr))/((V-b)^3)
    term_4 = (A4+B4*T+C4*exp(-k*Tr))/((V-b)^4)
    term_5 = (A5+B5*T+C5*exp(-k*Tr))/((V-b)^5)
    
    # cache -
    push!(term_array, term_1)
    push!(term_array, term_2)
    push!(term_array, term_3)
    push!(term_array, term_4)
    push!(term_array, term_5)

    # return -
    return sum(term_array)
end

function compute_temperature_mh(eos::MartinHouEquationOfState, P::Float64, V::Float64)::Float64 
    return 0.0
end
# ------------------------------------------------------------------------- #

# --- EXPORTED METHODS ----------------------------------------------------- #
function pressure(eos::AbstractEquationOfState, V::Float64, T::Float64)::Float64

    if typeof(eos) == IdealGasEquationOfState
        return compute_pressure_ideal_gas(eos,V,T)
    elseif (typeof(eos) == VanDerWaalsEquationOfState)
        return compute_pressure_vdw(eos,V,T)
    elseif (typeof(eos) == MartinHouEquationOfState)
        return compute_pressure_mh(eos, V, T)
    elseif (typeof(eos) == PengRobinsonEquationOfState)
        return compute_pressure_pr(eos, V, T)
    else
        throw(ArgumentError("$(typeof(eos)) is not supported"))
    end
end

function volume(eos::AbstractEquationOfState, P::Float64, T::Float64)::Float64
    
    if typeof(eos) == IdealGasEquationOfState
        return compute_volume_ideal_gas(eos, P, T)
    elseif (typeof(eos) == VanDerWaalsEquationOfState)
        return compute_volume_vdw(eos, P, T)
    elseif (typeof(eos) == MartinHouEquationOfState)
        return compute_volume_mh(eos, P, T)
    elseif (typeof(eos) == PengRobinsonEquationOfState)
        return compute_volume_pr(eos, P, T)
    else
        throw(ArgumentError("$(typeof(eos)) is not supported"))
    end
end

function temperature(eos::AbstractEquationOfState,P::Float64,V::Float64)::Float64
    if typeof(eos) == IdealGasEquationOfState
        return compute_temperature_ideal_gas(eos, P, V)
    elseif (typeof(eos) == VanDerWaalsEquationOfState)
        return compute_temperature_vdw(eos, P, V)
    elseif (typeof(eos) == MartinHouEquationOfState)
        return compute_temperature_mh(eos, P, V)
    elseif (typeof(eos) == PengRobinsonEquationOfState)
        return compute_temperature_pr(eos, P, V)
    else
        throw(ArgumentError("$(typeof(eos)) is not supported"))
    end
end
# ------------------------------------------------------------------------- #