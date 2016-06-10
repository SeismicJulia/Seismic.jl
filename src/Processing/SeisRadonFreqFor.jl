"""
    SeisRadonFreqFor(m; <keyword arguments>)

Transform a tau-p gather to a time-offset gather using a frequency domain
forward parabolic or linear Radon operator.

# Arguments
* `m::Array{T<:Real,2}`: 2D Radon panel, `m[1:nt,1:np]`, where `nt` is number of
time samples and `np` the number of curvatures or ray parameters.

# Keywords
* `order::ASCIIString="parab"`: `"parab"` for parabolic transform, `"linear"`
for linear transform.
* `dt::Real=0.002`: sampling interval in seconds.
* `h::Vector{Real}=collect(0.0:20.0:1000.0)`: offset vector.
* `href::Real=1000.0`: reference offset for parabolic Radon Transform.
* `p::Vector{Real}=collect(-0.05:0.01:2.2)`: if `order="parab"`, `p` is a vector
of residual moveout at reference offset in seconds ("curvatures");
if `order=linear`, `p` is a vector of ray parameters reference (s/m).
* `flow::Real=2.0`: minimum frequency in the data in Hz.
* `fhigh::Real=80.0`: maximum frequency in the data in Hz.

# Output
* `d`: data synthetized via forward Radon modeling, `d[1:nt, 1:nh]`.

# References
*  Hampson, D., 1986, Inverse velocity stacking for multiple elimination:
Canadian Journal of Exploration Geophysics, 22, 44-55.
*  Sacchi, M.D. and Ulrych, T.J., 1995, High-resolution velocity gathers and
offset space reconstruction: Geophysics, 60, 1169-1177.
"""
function SeisRadonFreqFor{Tm,Th,Tp}(m::Array{Tm,2}; order::ASCIIString="parab",
                                    dt::Real=0.002, href::Real=1000.0,
                                    h::Vector{Th}=collect(0.0:20.0:1000.0),
                                    p::Vector{Tp}=collect(-0.05:0.01:2.2),
                                    flow::Real=2.0, fhigh::Real=80.0)
    if order=="parab"
        I = 2
    elseif order=="linear"
        I = 1
        href = 1
    else
        error("Order should be equal to parab or linear")
    end
    nt, np = size(m)
    nw = nextpow2(nt)
    nh = length(h)
    m = cat(1, m, zeros(Tm, nw-nt, np))
    M = fft(m, 1)
    iw_low = round(Int, flow*dt*nw+1)
    iw_high = round(Int, fhigh*dt*nw+1)
    D = zeros(Complex64, nw, nh)
    for iw = iw_low:iw_high
        w  = 2.*pi*(iw-1)/(nw*dt)
        L = zeros(Complex64, nh, np)
        for ip = 1:np
            for ih = 1:nh
                phi = w*p[ip]*(h[ih]/href)^I
      	        L[ih, ip] = exp(-im*phi)
            end
        end 
        x = M[iw, :].'
        y = L*x
        D[iw, :] = y.'
    end 
    for iw = round(Int, nw/2)+2:nw
        D[iw, :] = conj(D[nw-iw+2, :])
    end
    d = real(ifft(D, 1))
    d = d[1:nt, :]
    return d
end
