module Seismic
    using Grid
    include("Utils/Utils.jl")
    include("Processing/Processing.jl")
    include("Imaging/Imaging.jl")
    include("Modelling/Modelling.jl")
    include("Solvers/Solvers.jl")
end
