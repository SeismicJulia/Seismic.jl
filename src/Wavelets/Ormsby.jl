"""
    Ormsby(; <keyword arguments>)

Create a Ormsby wavelet sampled every dt seconds with corner frequencies
defined by the vector f = [f1, f2, f3, f4]. The final wavelet is multiplied by
a Hamming window.

# Arguments

**Keyword arguments**

* `dt::Real=0.002`: sampling interval in secs.
* `f::Vector{Real}=[2.0, 10.0, 40.0, 60.0]`: corner frequencies in Hz.

      ^
    1 |     ***************
      |    *               *
      |   *                 *
      |  *                   *
      | *                     *
      -----------------------------> f
        f1  f2           f3  f4

# Example
```julia
julia> w = Ormsby(); plot(w);
```
"""

function Ormsby{T<:Real}(; dt::Real=0.002, f::Vector{T}=[2.0, 10.0, 40.0, 60.0])

    f1 = f[1]
    f2 = f[2]
    f3 = f[3]
    f4 = f[4]

    fc = (f2+f3)/2.0
    nw = 2.2/(fc*dt)
    nc = floor(Int, nw/2)
    t = dt*collect(-nc:1:nc)
    nw = 2*nc + 1
    a4 = (pi*f4)^2/(pi*(f4-f3))
    a3 = (pi*f3)^2/(pi*(f4-f3))
    a2 = (pi*f2)^2/(pi*(f2-f1))
    a1 = (pi*f1)^2/(pi*(f2-f1))

    u = a4*(sinc.(f4*t)).^2 - a3*(sinc.(f3*t)).^2
    v = a2*(sinc.(f2*t)).^2 - a1*(sinc.(f1*t)).^2

    w = u - v
    w = w.*Hamming(nw)/maximum(w)

end
