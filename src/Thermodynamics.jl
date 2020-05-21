module Thermodynamics

# include my codes (and external packages) -
include("./Include.jl")

# export methods -
export pressure
export temperature
export volume

# export types -
export AbstractEquationOfState
export VanDerWaalsEquationOfState
export IdealGasEquationOfState
export MartinHouEquationOfState

end # module
