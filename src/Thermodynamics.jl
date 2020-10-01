module Thermodynamics

# include my codes (and external packages) -
include("./Include.jl")

# export methods -
export pressure
export temperature
export volume
export compute_simple_compressibility_lk
export compute_complex_compressibility_lk

# factory methods -
export buildVanDerWaalsEquationOfState
export buildPengRobinsonEquationOfState
export buildLeeKeslerEquationOfState
export buildRedlichKwongEquationOfState
export buildSingleComponentWorkingFluid

# export types -
export AbstractEquationOfState
export AbstractWorkingFluid

export VanDerWaalsEquationOfState
export IdealGasEquationOfState
export MartinHouEquationOfState
export LeeKeslerEquationOfState
export RedlichKwongEquationOfState

export SingleComponentWorkingFluid
export MulticomponentWorkingFluid

end # module
