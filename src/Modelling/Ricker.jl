"""
    Ricker(; <keyword arguments>)

Ricker wavelet of central frequency f0 Hz sampled every dt seconds.

# Arguments

**Keyword arguments**

* `f0::Real=20.0`: central frequency in Hz.
* `dt::Real=0.002`: sampling interval in secs.

**Output**

* `w::Real`: Ricker wavelet.

# Examples
```julia
julia> w = Ricker(); plot(w);
```

**References**

* Sheriff, R., 2002, Encyclopedic dictionary of exploration geophysics, 4rd.
ed.: SEG. Geophysical Reference Series No. 13.

"""
function Ricker(; f0::Real=20.0, dt::Real=0.002)
    nw = 2.0/(f0*dt)
    nc = floor(Int, nw/2)
    t = dt*collect(-nc:1:nc)
    b = (pi*f0*t).^2
    w = (1.-2.*b).*exp(-b)
    return w
end
