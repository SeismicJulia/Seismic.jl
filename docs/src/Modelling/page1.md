# Generating Synthetic data set
Seismic.jl provides various ways to generate synthetic data set, like multi-dimensional linear, parabola, hyperbola events and finite-difference solver for acoustic wave equation (Currently only 2D is supported)


## SeisLinearEvents

```@docs
Seismic.SeisLinearEvents
```

## SeisAddNoise

```@docs
Seismic.SeisAddNoise
```

## SeisParabEvents

```@docs
Seismic.SeisParabEvents
```


## SeisHypEvents

```@docs
Seismic.SeisHypEvents
```


## SeisAddNoise

```@docs
Seismic.SeisAcousticWave 
```

### Example

```@example
using Seismic,PyPlot
SeisPlot(randn(20,20))
savefig("pic.svg"); nothing # hide
```
 ![](pic.svg)
