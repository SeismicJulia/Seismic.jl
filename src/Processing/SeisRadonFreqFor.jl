"""
    SeisRadonFreqFor(m, nt; <keyword arguments>)

Transform a tau-p gather to a time-offset gather using a frequency domain
forward parabolic or linear Radon operator.

# Arguments
* `m::Array{T<:Real,2}`: 2D Radon panel, `m[1:ntau,1:np]`, where `ntau` is the
number of intercept times and `np` the number of curvatures or ray parameters.
* `nt::Int`: number of time samples in the data domain.

# Keyword arguments
* `order::AbstractString="parab"`: `"parab"` for parabolic transform, `"linear"`
for linear transform.
* `dt::Real=0.004`: sampling interval in seconds.
* `h::Vector{Real}=collect(0.0:20.0:1000.0)`: offset vector; `h[1:nh]`.
* `href::Real=0.0`: reference offset for parabolic Radon Transform. If the
defautl value `href=0.0` is given, `href` is set to `max(abs(h))`.
* `p::Vector{Real}=collect(-0.05:0.01:2.2)`: `p[1:np]`. If `order="parab"`, `p`
is a vector of residual moveout ("curvatures") at reference offset `href` in
seconds; if `order=linear`, `p` is a vector of ray parameters in s/m.
* `flow::Real=0.0`: minimum frequency in the data in Hz.
* `fhigh::Real=125.0`: maximum frequency in the data in Hz.

# Output
* `d`: data synthetized via forward Radon modeling, `d[1:nt, 1:nh]`.

# References
*  Hampson, D., 1986, Inverse velocity stacking for multiple elimination:
Canadian Journal of Exploration Geophysics, 22, 44-55.
*  Sacchi, M.D. and Ulrych, T.J., 1995, High-resolution velocity gathers and
offset space reconstruction: Geophysics, 60, 1169-1177.
"""

function SeisRadonFreqFor{Tm<:Real, Th<:Real, Tp<:Real
                          }(m::Array{Tm,2}, nt::Int; order::AbstractString="parab",
                            dt::Real=0.004, href::Real=0.0,
                            h::Vector{Th}=collect(0.0:20.0:1000.0),
                            p::Vector{Tp}=collect(-0.05:0.01:2.2),
                            flow::Real=0.0, fhigh::Real=125.0)
    if order=="parab"
        I = 2
        href == 0 && (href = maximum(abs.(h)))
    elseif order=="linear"
        I = 1
        href = 1.0
    else
        error("Order should be equal to \"parab\" or \"linear\"")
    end
    ntau = size(m, 1)
    np = length(p)
    np == size(m, 2) || error("lenght(p) must be equal to size(m, 2)")
    nh = length(h)
    nw = 2*nextpow2(ntau)
    m = cat(1, m, zeros(Tm, nw-ntau, np))
    M = fft(m, 1)
    iw_low = round(Int, flow*dt*nw+1)
    iw_high = round(Int, fhigh*dt*nw+1)
    D = zeros(Complex{Tm}, nw, nh)
    for iw = iw_low:iw_high
        w  = 2.0*pi*(iw-1)/(nw*dt)
        L = zeros(Complex{Tm}, nh, np)
        for ip = 1:np
            for ih = 1:nh
                phi = w*p[ip]*(h[ih]/href)^I
      	        L[ih, ip] = exp(-im*phi)
            end
        end
        D[iw, :] = L*M[iw,:]
    end
    for iw = round(Int, nw/2)+2:nw
        D[iw, :] = conj(D[nw-iw+2, :])
    end
    d = real(ifft(D, 1))
    d = d[1:nt, :]
    return d
end
