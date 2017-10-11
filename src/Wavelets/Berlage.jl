"""
    Berlage(; <keyword arguments>)

Create a Berlage wavelet.

# Arguments

**Keyword arguments**

* `dt::Real=0.002`: sampling interval in secs.
* `f0::Real=20.0`: central frequency in Hz.
* `m::Real=2`: exponential parameter of Berlage wavelet.
* `alpha::Real=180.0`: alpha parameter of Berlage wavelet in rad/secs.
* `phi0::Real`: phase rotation in radians.

# Example
```julia
julia> w = Berlage(); plot(w);
```

**Reference**

* Aldridge, David F., 1990, The berlage wavelet: GEOPHYSICS, 55, 1508--1511.
"""

function Berlage(; dt::Real=0.002, f0::Real=20.0, m::Real=2, alpha::Real=180.0,
                 phi0::Real=0.0)

    nw = floor(Int, 2.2/(f0*dt))
    t = dt*collect(0:1:nw-1)
    w = (t.^m).*exp.(-alpha*t).*cos.(2*pi*f0*t + phi0);
    w = w/maximum(w)

end
