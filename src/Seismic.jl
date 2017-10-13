module Seismic
    using Interpolations,Requires,Compat
    include("ReadWrite/ReadWrite.jl")
    include("Utils/Utils.jl")
    include("Processing/Processing.jl")
    include("Imaging/Imaging.jl")
    include("Modelling/Modelling.jl")
    include("Operators/Operators.jl")
    include("Solvers/Solvers.jl")
    include("Tools/Tools.jl")
    include("Wavelets/Wavelets.jl")
    include("Windows/Windows.jl")
    @require PyPlot include("Plotting/Plotting.jl")
end
