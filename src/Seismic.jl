module Seismic
    using Grid,Requires,Compat
    include("Utils/Utils.jl")
    include("Processing/Processing.jl")
    include("Imaging/Imaging.jl")
    include("Modelling/Modelling.jl")
    include("Solvers/Solvers.jl")
    @require PyPlot include("Plotting/Plotting.jl")
end
