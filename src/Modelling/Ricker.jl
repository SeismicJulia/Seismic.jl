"""
    Ricker(; <keyword arguments>)

Create a Ricker wavelet.

# Arguments

**Keyword arguments**

* `dt::Real=0.002`: sampling interval in secs.
* `f0::Real=20.0`: central frequency in Hz.

# Example
```julia
julia> w = Ricker(); plot(w);
```

**Reference**

* Sheriff, Robert, 2002, Encyclopedic Dictionary of Applied Geophysics, fourth
ed.: Society of Exploration Geophysicists. Geophysical Reference Series No. 13.
"""

function Ricker(; dt::Real=0.002, f0::Real=20.0)

    nw = 2.0/(f0*dt)
    nc = floor(Int, nw/2)
    t = dt*collect(-nc:1:nc)
    b = (pi*f0*t).^2
    w = (1.-2.*b).*exp(-b)
    
end
