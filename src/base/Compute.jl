# --- PRIVATE METHODS ------------------------------------------------------ #
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

function compute_pressure_vdw(eos::VanDerWaalsEquationOfState, V::Float64, T::Float64)::Float64

    # check - is V = 0
    if (V==0.0)
        throw(ArgumentError("Volume == 0: attraction term is singular"))
    elseif ((V-b) == 0.0)
        throw(ArgumentError("Volume == b: singular repulsion term"))
    end

    # get R from the eos wrapper -
    R = eos.R
    a = eos.a
    b = eos.b

    # compute terms -
    repulsion = (R*T)/(V-b)
    attraction = (a^2)/(V^2)

    # return pressure -
    return (repulsion - attraction)
end

function compute_pressure_mh(eos::MartinHouEquationOfState, V::Float64, T::Float64)::Float64
    return 1.0
end
# ------------------------------------------------------------------------- #

# --- EXPORTED METHODS ----------------------------------------------------- #
function P(eos::AbstractEquationOfState, V::Float64, T::Float64)::Float64

    if typeof(eos) == IdealGasEquationOfState
        return compute_pressure_ideal_gas(eos,V,T)
    elseif (typeof(eos) == VanDerWaalsEquationOfState)
        return compute_pressure_vdw(eos,V,T)
    elseif (typeof(eos) == MartinHouEquationOfState)
        return compute_pressure_mh(eos, V, T)
    else
        throw(ArgumentError("$(typeof(eos)) is not supported"))
    end
end

function V(eos::AbstractEquationOfState,P::Float64,T::Float64)::Float64
    return 0.0
end

function T(eos::AbstractEquationOfState,P::Float64,V::Float64)::Float64
    return 0.0
end
# ------------------------------------------------------------------------- #