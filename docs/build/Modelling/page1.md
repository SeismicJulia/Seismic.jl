
<a id='Modelling-1'></a>

# Modelling


Seismic.jl provides various ways to generate synthetic data set, like multi-dimensional linear, parabola, hyperbola events and finite-difference solver for acoustic wave equation (Currently only 2D is supported)


<a id='SeisLinearEvents-1'></a>

## SeisLinearEvents

<a id='Seismic.SeisLinearEvents' href='#Seismic.SeisLinearEvents'>#</a>
**`Seismic.SeisLinearEvents`** &mdash; *Function*.



```
SeisLinearEvents(; <keyword arguments>)
```

Generate five dimensional data `d` consisting of linear events.

**Arguments**

**Keyword arguments**

  * `ot=0.0`: first sample for the time axis in secs.
  * `dt=0.004`: sampling interval in secs.
  * `nt=500`: number of time samples.
  * `ox1=0.0`: first sample for the first spatial dimension in meters.
  * `dx1=10.0`: sample interval for the first spatial dimension in meters.
  * `nx1=100`: number of samples for the first spatial dimension.
  * `ox2=0.0`: first sample for the second spatial dimension in meters.
  * `dx2=10.0`: sample interval for the second spatial dimension in meters.
  * `nx2=1`: number of samples for the second spatial dimension.
  * `ox3=0.0`: second sample for the third spatial dimension in meters.
  * `dx3=10.0`: sample interval for the third spatial dimension in meters.
  * `nx3=1`: number of samples for the third spatial dimension.
  * `ox4=0.0`: third sample for the fourth spatial dimension in meters.
  * `dx4=10.0`: sample interval for the fourth spatial dimension in meters.
  * `nx4=1`:number of samples for the fourth spatial dimension.
  * `tau=[1.0, 1.6]`: intercept traveltimes for each event.
  * `p1=[0.0000,-0.0001]`
  * `p2=[0.0003, 0.0002]`
  * `p3=[-0.0001,-0.0001]`
  * `p4=[0.0001,-0.0000]`
  * `amp=[1.0,-1.0]`: amplitudes for each linear event.
  * `f0=20.0`: central frequency of wavelet for each linear event.

**Example**

```julia
julia> d,extent = SeisLinearEvents(); SeisPlot(d);
```

Credits: Aaron Stanton, 2015


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Modelling/SeisLinearEvents.jl#L1-L39' class='documenter-source'>source</a><br>


<a id='SeisParabEvents-1'></a>

## SeisParabEvents

<a id='Seismic.SeisParabEvents' href='#Seismic.SeisParabEvents'>#</a>
**`Seismic.SeisParabEvents`** &mdash; *Function*.



```
SeisParabEvents(; <keyword arguments>)
```

Generate five dimensional data `d` consisting of parabolic events.

**Arguments**

**Keyword arguments**

  * `ot=0.0`: first sample for the time axis in secs.
  * `dt=0.004`: sampling interval in secs.
  * `nt=500`: number of time samples.
  * `ox1=0.0`: first sample for the first spatial dimension in meters.
  * `dx1=10.0`: sample interval for the first spatial dimension in meters.
  * `nx1=100`: number of samples for the first spatial dimension.
  * `ox2=0.0`: first sample for the second spatial dimension in meters.
  * `dx2=10.0`: sample interval for the second spatial dimension in meters.
  * `nx2=1`: number of samples for the second spatial dimension.
  * `ox3=0.0`: second sample for the third spatial dimension in meters.
  * `dx3=10.0`: sample interval for the third spatial dimension in meters.
  * `nx3=1`: number of samples for the third spatial dimension.
  * `ox4=0.0`: third sample for the fourth spatial dimension in meters.
  * `dx4=10.0`: sample interval for the fourth spatial dimension in meters.
  * `nx4=1`:number of samples for the fourth spatial dimension.
  * `tau=[1.0, 1.6]`: intercept traveltimes for each event.
  * `p1=[0.0000,-0.0001]
  * `p2=[0.0003, 0.0002]
  * `p3=[-0.0001,-0.0001]
  * `p4=[0.0001,-0.0000]
  * `amp=[1.0,-1.0]`: amplitudes for each parabolic event.
  * `wavelet="ricker"`: wavelet used to model the parabolicr events.
  * `f0=[20.0]`: central frequency of wavelet for each parabolic event.

**Example**

```julia
julia> d, extent = SeisParabEvents(); SeisPlot(d);
```

Credits: Mauricio D Sacchi, 2015


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Modelling/SeisParabEvents.jl#L1-L40' class='documenter-source'>source</a><br>


<a id='SeisHypEvents-1'></a>

## SeisHypEvents

<a id='Seismic.SeisHypEvents' href='#Seismic.SeisHypEvents'>#</a>
**`Seismic.SeisHypEvents`** &mdash; *Function*.



SeisHypEvents(; <keyword arguments>)

Generate two dimensional data `d` consisting of hyperbolic events.

**Keyword arguments**

  * `ot::Real=0.0`: first sample for the time axis in secs.
  * `dt::Real=0.004`: sampling interval in secs.
  * `nt::Int=301`: number of time samples.
  * `ox::Real=-1000.0`: first sample for spatial dimension in meters.
  * `dx::Real=20.0`: sample interval for the spatial dimension in meters.
  * `nx::Int=101`: number of samples for the spatial dimension.
  * `tau::Vector{Real}=[0.2, 0.6, 0.9]`: intercept traveltimes for each event.
  * `vel::Vector{Real}=[1500.0, 2000.0, 3000.0]`: rms velocities in m/s
  * `apex::Vector{Real}=[0.0, 0.0, 0.0]`: apex-shifts in meters.
  * `amp::Vector{Real}=[1.0, -1.0, 1.0]`: amplitudes for each event.
  * `wavelet::AbstractString="ricker"`: wavelet used to model the events.
  * `f0::Vector{Real}=[20.0]`: central frequency of wavelet for each event.

**Output**

  * `d::Array{Real, 2}`: two dimensional data consisting of hyperbolic events.
  * `extent::Extent`: extent of the data `d`.

**Examples**

```julia
julia> d, extent = SeisHypEvents(); SeisPlot(d, extent);
julia> d, extent = SeisHypEvents(apex=[100, 200, -300], f0=[30, 20, 15]);
SeisPlot(d, extent);
```


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Modelling/SeisHypEvents.jl#L1-L30' class='documenter-source'>source</a><br>


<a id='SeisAddNoise-1'></a>

## SeisAddNoise

<a id='Seismic.SeisAddNoise' href='#Seismic.SeisAddNoise'>#</a>
**`Seismic.SeisAddNoise`** &mdash; *Function*.



```
SeisAddNoise(d, snr; <keyword arguments>)
```

Add noise at a given signal-to-noise ratio level `snr` to an N-dimensional input data `d`. Noise can be band limited using kewyord `L`.

**Arguments**

  * `d::Array{Real, N}`: N-dimensional data.
  * `snr::Real`: signal-to-noise ratio.

**Keyword arguments**

  * `db::Bool=false`: `db=false` if `snr` is given by amplitude, `db=true` if

snr is given in dB.

  * `pdf::AbstractString="gaussian"`: random noise probability distribution:

`"gaussian"` or `"uniform"`.

  * `L::Int=1`: averaging operator length to band-limit the random noise.

**Examples**

```
julia> w = Ricker(); wn = SeisAddNoise(w, 2); plot(w); plot(wn);
MeasureSNR(w, wn)

julia> d, extent = SeisHypEvents(); dn = SeisAddNoise(d, 1.0, db=true, L=9);
SeisPlot([d dn], extent); MeasureSNR(d, dn, db=true)
```

Credits: Juan I. Sabbione, 2016


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Modelling/SeisAddNoise.jl#L1-L27' class='documenter-source'>source</a><br>


<a id='SeisAcousticWave-1'></a>

## SeisAcousticWave

<a id='Seismic.SeisAcousticWave' href='#Seismic.SeisAcousticWave'>#</a>
**`Seismic.SeisAcousticWave`** &mdash; *Function*.



```
shot = SeisAcousticWave(fidMtx, pos, isz, isx, f0, dt, tmax=2.0)
```

finite difference modeling of acoustic wave field, generate a common shot gather

**Arguments**

  * `fidMtx :: FidMtx`        : composite type of sparse matrix
  * `pos    :: Array{Int64,2}`: index of receivers, first column is the vertical index, second column contains horizontal index.
  * `isz    :: Int64`         : vertical index of source
  * `isz    :: Int64`         : horizontal index of source
  * `f0     :: Float64`       : dominant frequency of Ricker wavelet
  * `dt     :: Float64`       : size of time step

**keywords Arguments**

  * `tmax=1.0` : length of simulation

**Output**

  * `shot :: ShotGather`: composite type for common shot gather


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Modelling/SeisAcousticWave.jl#L445-L463' class='documenter-source'>source</a><br>

