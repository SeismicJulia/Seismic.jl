
<a id='Seismic.jl-1'></a>

# Seismic.jl


Seismic signal analysis and imaging


---


<a id='Overview-1'></a>

## Overview


Seismic.jl provides tools to **process**, **image**, and **plot** reflection seismic data in the Julia language.


<a id='_Convert-data-to-a-simple-format_-1'></a>

### _Convert data to a simple format_


Data and headers are stored separately as _filename.seisd_ and _filename.seish_ for simplicity. Functions are available to convert from popular formats such as SEGY, SU and RSF.


<a id='_Manipulate-data_-1'></a>

### _Manipulate data_


Contains functions for geometry calculation, sorting, windowing, patching/un-patching, and processing keyed on header word.


<a id='_Plot-publication-quality-figures_-1'></a>

### _Plot publication quality figures_


Produce color and wiggle plots using PyPlot.jl


---


<a id='Installation-1'></a>

## Installation


To install Seismic.jl type:


```no-highlight
Pkg.add("Seismic")
```


<a id='Index-1'></a>

## Index

- [`Seismic.Berlage`](Wavelets/page1.md#Seismic.Berlage)
- [`Seismic.Ormsby`](Wavelets/page1.md#Seismic.Ormsby)
- [`Seismic.Ricker`](Wavelets/page1.md#Seismic.Ricker)
- [`Seismic.SeisAcousticWave`](Modelling/page1.md#Seismic.SeisAcousticWave)
- [`Seismic.SeisAddNoise`](Modelling/page1.md#Seismic.SeisAddNoise)
- [`Seismic.SeisBinData`](Utils/page1.md#Seismic.SeisBinData)
- [`Seismic.SeisBinHeaders`](Utils/page1.md#Seismic.SeisBinHeaders)
- [`Seismic.SeisGeometry`](Utils/page1.md#Seismic.SeisGeometry)
- [`Seismic.SeisHypEvents`](Modelling/page1.md#Seismic.SeisHypEvents)
- [`Seismic.SeisLinearEvents`](Modelling/page1.md#Seismic.SeisLinearEvents)
- [`Seismic.SeisParabEvents`](Modelling/page1.md#Seismic.SeisParabEvents)
- [`Seismic.SeisPatch`](Utils/page1.md#Seismic.SeisPatch)
- [`Seismic.SeisPlot`](Plotting/page1.md#Seismic.SeisPlot)
- [`Seismic.SeisPlotCoordinates`](Plotting/page1.md#Seismic.SeisPlotCoordinates)
- [`Seismic.SeisSort`](Utils/page1.md#Seismic.SeisSort)
- [`Seismic.SeisUnPatch`](Utils/page1.md#Seismic.SeisUnPatch)
- [`Seismic.SeisWindow`](Utils/page1.md#Seismic.SeisWindow)

