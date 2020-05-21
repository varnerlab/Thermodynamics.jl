# base equation of state type -
abstract type AbstractEquationOfState end

mutable struct VanDerWaalsEquationOfState <: AbstractEquationOfState

    R::Float64
    a::Float64
    b::Float64

    function VanDerWaalsEquationOfState(R::Float64, a::Float64, b::Float64)
        this = new(R,a,b)
    end
end

mutable struct IdealGasEquationOfState <: AbstractEquationOfState

    # parameters -
    R::Float64

    function IdealGasEquationOfState(R::Float64)
        this = new(R)
    end
end

mutable struct MartinHouEquationOfState <: AbstractEquationOfState

    # parameters -
    parameters::Dict{String,Any}

    function MartinHouEquationOfState(parameterDict::Dict{String, Any})
        this = new(parameterDict)
    end
end
