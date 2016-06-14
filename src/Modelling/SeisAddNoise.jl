"""
    SeisAddNoise(d, snr; <keyword arguments>)

Add band-limited gaussian noise at a given signal-to-noise ratio level `snr` to
a N-dimensional input data `d`.

# Arguments
* `d::Array{Real, N}`: N-dimensional data.
* `snr::Real`: signal-to-noise ratio

# Keyword arguments
* `noise::ASCIIString="gaussian"`: random noise distribution. 
* `L::Int=9`: averaging operator length for random noise.

# Output
* `dnoisy::Array{Real, 2}`: data contaminated with noise.

# Example
```julia
julia> d, extent = SeisHypEvents(); dnoisy = SeisAddNoise(d, 2.0);
SeisPlot(dnoisy, extent);

```
"""

function SeisAddNoise{T<:Real, N}(d::Array(T, N}, snr::Real;
                                  , dt::Real=0.004, nt::Int=301,

