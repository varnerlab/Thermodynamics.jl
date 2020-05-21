module Thermodynamics

# include my codes (and external packages) -
include("./Include.jl")

# export methods -
export P
export T
export V

# export types -
export AbstractEquationOfState
export VanDerWaalsEquationOfState
export IdealGasEquationOfState
export MartinHouEquationOfState

end # module
