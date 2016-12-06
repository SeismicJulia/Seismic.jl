"""
    MeasureSNR(signal, noisy; db=false)

Measure the signal-to-noise ratio between the clean input `signal` and the
contaminated input `noisy`.

# Arguments
* `signal::Array{Real, N}`: N-dimensional clean signal. `N` must be <= 5.
* `noisy::Array{Real, N}`: N-dimensional noisy signal of same size as `signal`.

# Keyword arguments
* `db::Bool=false`: `db=false` if the signal-to-noise ratio is measured by
amplitude, or `db=true` if the signal-to-noise ratio is measure in dB.

# Example
```julia
julia> d, extent = SeisHypEvents(); dnoisy = SeisAddNoise(d, 2); 
MeasureSNR(d, dnoisy)
```
"""

function MeasureSNR{Ts<:Real, Tn<:Real, N}(signal::Array{Ts, N},
                                           noisy::Array{Tn, N}; db::Bool=false)

    if db==false
        snr = vecnorm(signal)/vecnorm(signal-noisy)
    elseif db==true
        snr = 20*log10(vecnorm(signal)/vecnorm(signal-noisy))
    end
    return snr

end
