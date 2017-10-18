
<a id='Wavelets-1'></a>

# Wavelets


A collection of seismic data processing wavelet


<a id='Berlage-1'></a>

## Berlage

<a id='Seismic.Berlage' href='#Seismic.Berlage'>#</a>
**`Seismic.Berlage`** &mdash; *Function*.



```
Berlage(; <keyword arguments>)
```

Create a Berlage wavelet.

**Arguments**

**Keyword arguments**

  * `dt::Real=0.002`: sampling interval in secs.
  * `f0::Real=20.0`: central frequency in Hz.
  * `m::Real=2`: exponential parameter of Berlage wavelet.
  * `alpha::Real=180.0`: alpha parameter of Berlage wavelet in rad/secs.
  * `phi0::Real`: phase rotation in radians.

**Example**

```julia
julia> w = Berlage(); plot(w);
```

**Reference**

  * Aldridge, David F., 1990, The berlage wavelet: GEOPHYSICS, 55, 1508–1511.


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Wavelets/Berlage.jl#L1-L24' class='documenter-source'>source</a><br>


<a id='Ormsby-1'></a>

## Ormsby

<a id='Seismic.Ormsby' href='#Seismic.Ormsby'>#</a>
**`Seismic.Ormsby`** &mdash; *Function*.



```
Ormsby(; <keyword arguments>)
```

Create a Ormsby wavelet sampled every dt seconds with corner frequencies defined by the vector f = [f1, f2, f3, f4]. The final wavelet is multiplied by a Hamming window.

**Arguments**

**Keyword arguments**

  * `dt::Real=0.002`: sampling interval in secs.
  * `f::Vector{Real}=[2.0, 10.0, 40.0, 60.0]`: corner frequencies in Hz.

    ```
    ^
    ```

    1 |     ***************     |    *               *     |   *                 *     |  *                   *     | *                     *     ––––––––––––––-> f       f1  f2           f3  f4

**Example**

```julia
julia> w = Ormsby(); plot(w);
```


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Wavelets/Ormsby.jl#L1-L28' class='documenter-source'>source</a><br>


<a id='Ricker-1'></a>

## Ricker

<a id='Seismic.Ricker' href='#Seismic.Ricker'>#</a>
**`Seismic.Ricker`** &mdash; *Function*.



```
Ricker(; <keyword arguments>)
```

Create a Ricker wavelet.

**Keyword arguments**

  * `dt::Real=0.002`: sampling interval in secs.
  * `f0::Real=20.0`: central frequency in Hz.

**Examples**

```julia
julia> w = Ricker(); plot(w);
julia> w = Ricker(dt=0.004, f0=20); plot(w);
```

**Reference**

Sheriff, Robert, 2002, Encyclopedic Dictionary of Applied Geophysics, fourth ed.: Society of Exploration Geophysicists. Geophysical Reference Series No. 13.


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Wavelets/Ricker.jl#L1-L19' class='documenter-source'>source</a><br>

