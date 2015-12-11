using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.lines as lines
export SeisPlot,
SeisPlotFKSpectrum,
SeisPlotAmplitudeSpectrum,
SeisMap
include("SeisPlot.jl")
include("SeisPlotFKSpectrum.jl")
include("SeisPlotAmplitudeSpectrum.jl")
include("SeisMap.jl")

