# base equation of state type -
abstract type AbstractEquationOfState end
abstract type AbstractWorkingFluid end

mutable struct SingleComponentWorkingFluid <: AbstractWorkingFluid

    # data -
    Tc::Float64 # units: K
    Pc::Float64 # units: MPa
    Vc::Float64 # units: L
    ω::Float64  # units: dimensionless
    R::Float64  # units: MPa L mol^-1 K^-1 

    function SingleComponentWorkingFluid(Tc::Float64, Pc::Float64, Vc::Float64, ω::Float64; R=0.008314)
        this = new(Tc,Pc,Vc,ω,R)
    end
end

mutable struct MulticomponentWorkingFluid <: AbstractWorkingFluid
end

mutable struct RedlichKwongEquationOfState <: AbstractEquationOfState

    R::Float64
    a::Float64
    b::Float64

    function RedlichKwongEquationOfState(R::Float64, a::Float64, b::Float64)
        this = new(R,a,b)
    end
end

mutable struct VanDerWaalsEquationOfState <: AbstractEquationOfState

    R::Float64
    a::Float64
    b::Float64

    function VanDerWaalsEquationOfState(R::Float64, a::Float64, b::Float64)
        this = new(R,a,b)
    end
end

mutable struct LeeKeslerEquationOfState <: AbstractEquationOfState

    # parameters -
    parameters::Dict{String,Any}

    function LeeKeslerEquationOfState(parameterDict::Dict{String, Any})
        this = new(parameterDict)
    end
end

mutable struct PengRobinsonEquationOfState <: AbstractEquationOfState

    # parameters -
    R::Float64
    a::Float64
    b::Float64
    𝞳::Float64
    Tc::Float64

    function PengRobinsonEquationOfState(R::Float64,a::Float64, b::Float64, 𝞳::Float64, Tc::Float64)
        this = new(R,a,b,𝞳,Tc)
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
