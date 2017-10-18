
<a id='Tools-1'></a>

# Tools


<a id='MeasureSNR-1'></a>

## MeasureSNR

<a id='Seismic.MeasureSNR' href='#Seismic.MeasureSNR'>#</a>
**`Seismic.MeasureSNR`** &mdash; *Function*.



```
MeasureSNR(signal, noisy; db=false)
```

Measure the signal-to-noise ratio between the clean input `signal` and the contaminated input `noisy`.

**Arguments**

  * `signal::Array{Real, N}`: N-dimensional clean signal. `N` must be <= 5.
  * `noisy::Array{Real, N}`: N-dimensional noisy signal of same size as `signal`.

**Keyword arguments**

  * `db::Bool=false`: `db=false` if the signal-to-noise ratio is measured by

amplitude, or `db=true` if the signal-to-noise ratio is measure in dB.

**Example**

```julia
julia> d, extent = SeisHypEvents(); dnoisy = SeisAddNoise(d, 2); 
MeasureSNR(d, dnoisy)
```


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Tools/MeasureSNR.jl#L1-L20' class='documenter-source'>source</a><br>


<a id='PadFirstAxis-1'></a>

## PadFirstAxis

<a id='Seismic.PadFirstAxis' href='#Seismic.PadFirstAxis'>#</a>
**`Seismic.PadFirstAxis`** &mdash; *Function*.



```
PadFirstAxis(d, N1)
```

Zero-padding first axis of N-dimensional array `d`.

**Arguments**

  * `d::Array{Real, N}`: N-dimensional data.
  * `N1::Int`: total number of samples for the first axis of the zero-padded data.

**Examples**

```julia
julia> w = Ricker(); wpad = PadFirstAxis(w, 128); plot(w); plot(wpad+1.0)

julia> d, ext = SeisHypEvents(); dpad = PadFirstAxis(d, 512);
SeisPlot(d); SeisPlot(dpad)
```


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Tools/PadFirstAxis.jl#L1-L17' class='documenter-source'>source</a><br>

