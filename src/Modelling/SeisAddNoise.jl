"""
    SeisAddNoise(d, snr; <keyword arguments>)

Add band-limited noise at a given signal-to-noise ratio level `snr` to
a N-dimensional input data `d`.

# Arguments
* `d::Array{Real, N}`: N-dimensional data.
* `snr::Real`: signal-to-noise ratio

# Keyword arguments
* `noise::ASCIIString="gaussian"`: random noise distribution: `"gaussian"` or
`"uniform"`. 
* `L::Int=9`: averaging operator length for random noise.
* `db::Bool=false`: `db=false` for `snr` is given by amplitude, `db=false` if
snr is given in dB.

# Examples
```julia
julia> w = Ricker(); wnoisy = SeisAddNoise(w, 2.0); plot(w); plot(wnoisy);
julia> d, extent = SeisHypEvents(); dnoisy = SeisAddNoise(d, 2.0);
SeisPlot([d dnoisy], extent);
```
"""

function SeisAddNoise{T<:Real, N}(d::Array(T, N}, snr::Real; 
                                  dt::Real=0.004, nt::Int=301,

