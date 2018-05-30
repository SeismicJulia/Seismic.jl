"""
```
SeisRadonFreqInv(d; <keyword arguments>)
```

Transform a Gather from time-offset gather to tau-p gather using a frequency
domain inverse parabolic or linear Radon operator via least-squares inversion.

# Arguments
* `d::Array{T<:Real,2}`: 2D data, `d[1:nt,1:nh]`, where `nt` is number of
time samples and `nh` the number of receivers.

# Keyword arguments
* `order::AbstractString="parab"`: `"parab"` for parabolic transform, `"linear"`
for linear transform.
* `dt::Real=0.004`: sampling interval in seconds.
* `h::Vector{Real}=collect(0.0:20.0:1000.0)`: offset vector `h[1:nh]`.
* `href::Real=0.0`: reference offset for parabolic Radon Transform. If the
defautl value `href=0.0` is given, `href` is set to `max(abs(h))`.
* `p::Vector{Real}=collect(-0.05:0.01:2.2)`: `p[1:np]`. If `order="parab"`, `p`
is a vector of residual moveout ("curvatures") at reference offset `href` in
seconds; if `order=linear`, `p` is a vector of ray parameters in s/m.
* `flow::Real=0.0`: minimum frequency in the data in Hz.
* `fhigh::Real=125.0`: maximum frequency in the data in Hz.
* `mu::Real=0.001`: trade off parameter or damping for the L.S. inversion.

# Output
* `m`: inverted Radon panel `m[1:ntau, 1:np]`.

# References
*  Hampson, D., 1986, Inverse velocity stacking for multiple elimination:
Canadian Journal of Exploration Geophysics, 22, 44-55.
*  Sacchi, M.D. and Ulrych, T.J., 1995, High-resolution velocity gathers and
offset space reconstruction: Geophysics, 60, 1169-1177.
"""

function SeisRadonFreqInv{Td<:Real, Th<:Real, Tp<:Real
                          }(d::Array{Td,2}; order::AbstractString="parab",
                            dt::Real=0.004, href::Real=0.0,
                            h::Vector{Th}=collect(0.0:20.0:1000.0),
                            p::Vector{Tp}=collect(-0.05:0.01:2.2),
                            flow::Real=0.0, fhigh::Real=125.0, mu::Real=0.001)
    if order=="parab"
        I = 2
        href == 0 && (href = maximum(abs.(h)))
    elseif order=="linear"
        I = 1
        href = 1.0
    else
        error("Order should be equal to \"parab\" or \"linear\"")
    end
    nt = size(d, 1)
    nh = length(h)
    nh == size(d, 2) || error("lenght(h) must be equal to size(d, 2)")
    np = length(p)
    nw = 2*nextpow2(nt)
    d = cat(1, d, zeros(Td, nw-nt, nh))
    D = fft(d, 1)
    iw_low = round(Int, flow*dt*nw+1)
    iw_high = round(Int, fhigh*dt*nw+1)
    M = zeros(Complex{Td}, nw, np)
    for iw = iw_low:iw_high
        w  = 2.0*pi*(iw-1)/(nw*dt)
        L = zeros(Complex{Td}, nh, np)
        for ip = 1:np
            for ih = 1:nh
                phi = w*p[ip]*(h[ih]/href)^I
      	        L[ih, ip] = exp(-im*phi)
            end
        end
        R = L'*L
        r0 = trace(R)/np
        xa = L'*D[iw, :]
        Q = mu*r0*eye(np)
        x  = (R + Q)\xa
        M[iw, :] = x.'
    end
    for iw = round(Int, nw/2)+2:nw
        M[iw, :] = conj(M[nw-iw+2, :])
    end
    m = real(ifft(M, 1))
    m = m[1:nt, :]
    return m
end
