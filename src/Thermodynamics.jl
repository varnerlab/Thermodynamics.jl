module Thermodynamics

# include my codes (and external packages) -
include("./Include.jl")

# export methods -
export pressure
export temperature
export volume

# factory methods -
export buildVanDerWaalsEquationOfState

# export types -
export AbstractEquationOfState
export AbstractWorkingFluid

export VanDerWaalsEquationOfState
export IdealGasEquationOfState
export MartinHouEquationOfState

export SingleComponentWorkingFluid
export MulticomponentWorkingFluid

end # module
