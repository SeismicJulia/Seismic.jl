
<a id='Processing-1'></a>

# Processing


<a id='SeisBandPass-1'></a>

## SeisBandPass

<a id='Seismic.SeisBandPass' href='#Seismic.SeisBandPass'>#</a>
**`Seismic.SeisBandPass`** &mdash; *Function*.



```
SeisBandPass(d ; <keyword arguments>)
```

Apply a bandpass filter to a 2D array input. Input and output data is in tx domain. Filter is applied in fx domain.

**Arguments**

  * `d`: Input 2D data array in tx domain.

**Keyword arguments**

  * `dt=0.001`: time sampling interval.
  * `fa=0,fb=0,fc=60,fd=80`: corner frequencies in Hz.

**Output**

  * `d`: Filtered 2d data arrayn in tx domain.

**Example**

```julia julia> d = SeisLinearEvents(d); SeisPlot(d,plot_type="Amplitude",dy=0.004,fmax=125); julia> d_filter = SeisBandPass(d;dt=0.004,fc=15,fd=35); SeisPlot(d_filter,plot_type="Amplitude",dy=0.004,fmax=125);


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisBandPass.jl#L1-L21' class='documenter-source'>source</a><br>


<a id='SeisDecimate-1'></a>

## SeisDecimate

<a id='Seismic.SeisDecimate' href='#Seismic.SeisDecimate'>#</a>
**`Seismic.SeisDecimate`** &mdash; *Function*.



```
SeisDecimate(d ; <keyword arguments>)
```

Decimate a multidimensional array input. Input and output have the same dimension

**Arguments**

  * `d`: Input data. Between 2D and 5D

**Keyword arguments**

  * `mode=random`: Data is decimated randomly. Otherwise it is regularly decimated.
  * `perc=50`: percentage of traces to be decimated from original array
  * `incx1=1,incx2=1,incx3=1,incx4=1`: traces to decimate in case of regular decimation

**Output**

  * `out`: Decimated data.

**Example**

```julia julia> d = SeisLinearEvents(d); SeisPlot(d) julia> d_dec = SeisDecimate(d); SeisPlot(d_dec)


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisDecimate.jl#L1-L21' class='documenter-source'>source</a><br>


<a id='SeisDelay-1'></a>

## SeisDelay

<a id='Seismic.SeisDelay' href='#Seismic.SeisDelay'>#</a>
**`Seismic.SeisDelay`** &mdash; *Function*.



```
SeisDelay(d ; <keyword arguments>)
```

Apply a time delay to a 2D array input. 

**Arguments**

  * `d`: Input 2D data array in tx domain.

**Keyword arguments**

  * `delay=0.1`: time delay in seconds.
  * `dt=0.001`: time sampling

**Output**

  * `d2`: Delayed data in time domain

**Example**

```julia julia> d = SeisLinearEvents(d); SeisPlot(d); julia> d2 = SeisDelay(d;dt=0.004); SeisPlot(d);


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisDelay.jl#L1-L21' class='documenter-source'>source</a><br>


<a id='SeisEnvelope-1'></a>

## SeisEnvelope

<a id='Seismic.SeisEnvelope' href='#Seismic.SeisEnvelope'>#</a>
**`Seismic.SeisEnvelope`** &mdash; *Function*.



```
SeisEnvelope(d)
```

Calculate the envelope attribute of an input trace

**Arguments**

  * `d`: Input data.

**Output**

  * `out`: Envelope of input data

**Example**

```julia julia> d = SeisLinearEvents(d); SeisPlot(d) julia> out = SeisEnvelope(d); SeisPlot(d_dec)


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisEnvelope.jl#L1-L16' class='documenter-source'>source</a><br>


<a id='SeisFKFilter-1'></a>

## SeisFKFilter

<a id='Seismic.SeisFKFilter' href='#Seismic.SeisFKFilter'>#</a>
**`Seismic.SeisFKFilter`** &mdash; *Function*.



```
SeisFKFilter(d ; <keyword arguments>)
```

Decimate a multidimensional array input. Input and output have the same dimension

**Arguments**

  * `d`: 2D Input data.

**Keyword arguments**

  * `dt=0.002`: time sampling interval.
  * `dx=10`: space sampling interval.
  * `va=-2000,vb=-3000,vc=3000,vd=2000`: corner velocities to be filtered

**Output**

  * `out`: Filtered data.


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisFKFilter.jl#L1-L18' class='documenter-source'>source</a><br>


<a id='SeisGain-1'></a>

## SeisGain

<a id='Seismic.SeisGain' href='#Seismic.SeisGain'>#</a>
**`Seismic.SeisGain`** &mdash; *Function*.



```
SeisGain(d, t; <keyword arguments>)
```

Gain a group of traces.

**Arguments**

  * `d::Array{Real,2}`: two dimensional data.

**Keyword arguments**

  * `dt::Real=0.002`: sampling interval in secs.
  * `kind::AbstractString="time"`: if kind="time", gain = t.^a . * exp(-bt);                        if kind="agc", automatic gain control is applied.
  * `param::Vector{Real}=[2.0,0.0]`: if kind="time", param = [a,b];                                  if kind="agc", param = [agc_gate]
  * `norm::Int=0`: `norm=0` no normalization; `norm=1` normalize each trace by                 amplitude; `norm=2` normalize each trace by rms value/

**Output**

  * d::Array{Real, 2}`: gained two dimensional data.

**Example**

```julia
julia> d, extent = SeisHypEvents();
       dout = SeisGain(d, kind="agc", param=[0.05]);
       SeisPlot([d dout], extent);
```

Credits: Juan I. Sabbione, Aaron Staton, Mauricio D. Sacchi, 2016


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisGain.jl#L1-L30' class='documenter-source'>source</a><br>


<a id='SeisKolmogoroff-1'></a>

## SeisKolmogoroff

<a id='Seismic.SeisKolmogoroff' href='#Seismic.SeisKolmogoroff'>#</a>
**`Seismic.SeisKolmogoroff`** &mdash; *Function*.



```
SeisKolmogoroff(w)
```

Kolmogoroff factorization. Transform a wavelet into its minimum phase equivalent.

**Arguments**

  * `w::Real`: input wavelet.

**Example**

```julia
julia> w = Ricker()
julia> wmin = SeisKolmogoroff(w)
julia> plot(w); plot(wmin)
```

**Reference**

  * Claerbout, Jon F., 1976, Fundamentals of geophysical data processing.

McGraw-Hill Inc.


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisKolmogoroff.jl#L1-L21' class='documenter-source'>source</a><br>


<a id='SeisMWNI-1'></a>

## SeisMWNI

<a id='Seismic.SeisMWNI' href='#Seismic.SeisMWNI'>#</a>
**`Seismic.SeisMWNI`** &mdash; *Function*.



**SeisMWNI**

*Minimum Weighted Norm Interpolation of seismic records.*

**IN**   

  * d_in: input data that can have up to 5 dimensions
  * dt=0.001 sampling rate along the time axis (in seconds)
  * fmax=99999. maximum temporal frequency to process.
  * padt=2 padding to use for the time axis
  * padx=1 padding to use for the spatial axes
  * Niter_internal=10 number of internal iterations for Conjugate Gradients
  * Niter_external=3 number of external iterations for iterative reweighting

**OUT**  

  * d_out: interpolated data


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisMWNI.jl#L1-L19' class='documenter-source'>source</a><br>


<a id='SeisPOCS-1'></a>

## SeisPOCS

<a id='Seismic.SeisPOCS' href='#Seismic.SeisPOCS'>#</a>
**`Seismic.SeisPOCS`** &mdash; *Function*.



**SeisPOCS**

*Projection Onto Convex Sets interpolation of seismic records.*

**IN**

  * d_in: input data that can have up to 5 dimensions
  * p=1., exponent for thresholding (1 is equivalent to soft thres. high number is equivalent to hard thresholding)
  * alpha=1 add-back ratio for imputation step. Use 1 for noise free data, and < 1 for denoising of original traces.
  * dt=0.001 sampling rate along the time axis (in seconds)
  * fmax=99999. maximum temporal frequency to process.
  * padt=2 padding to use for the time axis
  * padx=1 padding to use for the spatial axes
  * Niter=100 number of iterations

**OUT**

  * d_out: interpolated data


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisPOCS.jl#L1-L20' class='documenter-source'>source</a><br>


<a id='SeisRadonFreqFor-1'></a>

## SeisRadonFreqFor

<a id='Seismic.SeisRadonFreqFor' href='#Seismic.SeisRadonFreqFor'>#</a>
**`Seismic.SeisRadonFreqFor`** &mdash; *Function*.



```
SeisRadonFreqFor(m, nt; <keyword arguments>)
```

Transform a tau-p gather to a time-offset gather using a frequency domain forward parabolic or linear Radon operator.

**Arguments**

  * `m::Array{T<:Real,2}`: 2D Radon panel, `m[1:ntau,1:np]`, where `ntau` is the

number of intercept times and `np` the number of curvatures or ray parameters.

  * `nt::Int`: number of time samples in the data domain.

**Keyword arguments**

  * `order::AbstractString="parab"`: `"parab"` for parabolic transform, `"linear"`

for linear transform.

  * `dt::Real=0.004`: sampling interval in seconds.
  * `h::Vector{Real}=collect(0.0:20.0:1000.0)`: offset vector; `h[1:nh]`.
  * `href::Real=0.0`: reference offset for parabolic Radon Transform. If the

defautl value `href=0.0` is given, `href` is set to `max(abs(h))`.

  * `p::Vector{Real}=collect(-0.05:0.01:2.2)`: `p[1:np]`. If `order="parab"`, `p`

is a vector of residual moveout ("curvatures") at reference offset `href` in seconds; if `order=linear`, `p` is a vector of ray parameters in s/m.

  * `flow::Real=0.0`: minimum frequency in the data in Hz.
  * `fhigh::Real=125.0`: maximum frequency in the data in Hz.

**Output**

  * `d`: data synthetized via forward Radon modeling, `d[1:nt, 1:nh]`.

**References**

  * Hampson, D., 1986, Inverse velocity stacking for multiple elimination:

Canadian Journal of Exploration Geophysics, 22, 44-55.

  * Sacchi, M.D. and Ulrych, T.J., 1995, High-resolution velocity gathers and

offset space reconstruction: Geophysics, 60, 1169-1177.


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisRadonFreqFor.jl#L1-L33' class='documenter-source'>source</a><br>


<a id='SeisRadonFreqInv-1'></a>

## SeisRadonFreqInv

<a id='Seismic.SeisRadonFreqInv' href='#Seismic.SeisRadonFreqInv'>#</a>
**`Seismic.SeisRadonFreqInv`** &mdash; *Function*.



```
SeisRadonFreqInv(d; <keyword arguments>)
```

Transform a CMP Gather time-offset gather to tau-p gather using a frequency domain inverse parabolic or linear Radon operator via least-squares inversion.

**Arguments**

  * `d::Array{T<:Real,2}`: 2D data, `d[1:nt,1:nh]`, where `nt` is number of

time samples and `nh` the number of receivers.

**Keyword arguments**

  * `order::AbstractString="parab"`: `"parab"` for parabolic transform, `"linear"`

for linear transform.

  * `dt::Real=0.004`: sampling interval in seconds.
  * `h::Vector{Real}=collect(0.0:20.0:1000.0)`: offset vector `h[1:nh]`.
  * `href::Real=0.0`: reference offset for parabolic Radon Transform. If the

defautl value `href=0.0` is given, `href` is set to `max(abs(h))`.

  * `p::Vector{Real}=collect(-0.05:0.01:2.2)`: `p[1:np]`. If `order="parab"`, `p`

is a vector of residual moveout ("curvatures") at reference offset `href` in seconds; if `order=linear`, `p` is a vector of ray parameters in s/m.

  * `flow::Real=0.0`: minimum frequency in the data in Hz.
  * `fhigh::Real=125.0`: maximum frequency in the data in Hz.
  * `mu::Real=0.001`: trade off parameter or damping for the L.S. inversion.

**Output**

  * `m`: inverted Radon panel `m[1:ntau, 1:np]`.

**References**

  * Hampson, D., 1986, Inverse velocity stacking for multiple elimination:

Canadian Journal of Exploration Geophysics, 22, 44-55.

  * Sacchi, M.D. and Ulrych, T.J., 1995, High-resolution velocity gathers and

offset space reconstruction: Geophysics, 60, 1169-1177.


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Processing/SeisRadonFreqInv.jl#L1-L33' class='documenter-source'>source</a><br>

