using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.lines as lines
export SeisPlot,
SeisPlotFKSpectrum,
SeisPlotAmplitudeSpectrum,
SeisPlotCoordinates
include("SeisPlot.jl")
include("SeisPlotFKSpectrum.jl")
include("SeisPlotAmplitudeSpectrum.jl")
include("SeisPlotCoordinates.jl")

