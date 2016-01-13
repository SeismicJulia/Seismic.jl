module Seismic
    using Lexicon,Docile,Grid,Requires,Compat
    @compat include("Utils/Utils.jl")
    @compat include("Processing/Processing.jl")
    @compat include("Imaging/Imaging.jl")
    @compat include("Modelling/Modelling.jl")
    @compat include("Solvers/Solvers.jl")
    @require PyPlot include("Plotting/Plotting.jl")
end
