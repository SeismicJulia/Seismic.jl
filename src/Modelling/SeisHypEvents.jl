"""
   SeisHypEvents(; <keyword arguments>)

Generate two dimensional data `d` consisting of hyperbolic events.

# Keyword arguments
* `ot::Real=0.0`: first sample for the time axis in secs.
* `dt::Real=0.004`: sampling interval in secs.
* `nt::Int=301`: number of time samples.
* `ox::Real=-1000.0`: first sample for spatial dimension in meters.
* `dx::Real=20.0`: sample interval for the spatial dimension in meters.
* `nx::Int=101`: number of samples for the spatial dimension.
* `tau::Vector{Real}=[0.2, 0.6, 0.9]`: intercept traveltimes for each event.
* `vel::Vector{Real}=[1500.0, 2000.0, 3000.0]`: rms velocities in m/s
* `apex::Vector{Real}=[0.0, 0.0, 0.0]`: apex-shifts in meters.
* `amp::Vector{Real}=[1.0, -1.0, 1.0]`: amplitudes for each event.
* `wavelet::AbstractString="ricker"`: wavelet used to model the events.
* `f0::Vector{Real}=[20.0]`: central frequency of wavelet for each event.

# Output
* `d::Array{Real, 2}`: two dimensional data consisting of hyperbolic events.
* `extent::Extent`: extent of the data `d`.

# Examples
```julia
julia> d, extent = SeisHypEvents(); SeisPlot(d, extent);
julia> d, extent = SeisHypEvents(apex=[100, 200, -300], f0=[30, 20, 15]);
SeisPlot(d, extent);
```
"""

function SeisHypEvents{T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real
                      }(; ot::Real=0.0, dt::Real=0.004, nt::Int=301,
                        ox::Real=-1000.0, dx::Real=20.0, nx::Int=101,
                        tau::Vector{T1}=[0.2, 0.6, 0.9],
                        vel::Vector{T2}=[1500.0, 2000.0, 3000.0],
                        apex::Vector{T3}=[0.0, 0.0, 0.0],
                        amp::Vector{T4}=[1.0, -1.0, 1.0],
                        wavelet::AbstractString="ricker",
                        f0::Vector{T5}=[20.0])

    x = ox + collect(0:1:nx-1)*dx
    if length(f0) != length(tau)
	f0 = f0[1] + 0.0*tau
    end
    nf = 4*nextpow2(nt)
    dw = 2.0*pi/(nf*dt)
    nw = floor(Int, floor(nf/2)) + 1
    D = zeros(Complex{Float64}, nf, nx)
    nevents = length(tau)
    for ievent = 1:nevents
        if wavelet=="ricker"
            wav = Ricker(f0=f0[ievent], dt=dt)
            nwav = length(wav)
            wav = cat(1, wav, zeros(nf-nwav))
            Wav = fft(wav)
            delay = dt*(floor(Int, nwav/2))
            shift = sqrt.(tau[ievent]^2 + ((x-apex[ievent])/vel[ievent]).^2)'
        else
            error("Other wavelets apart from Ricker not implemented yet")
        end
        for iw = 1:nw
            w = (iw - 1)*dw
            D[iw:iw, :] += amp[ievent]*Wav[iw]*exp.(-1im*w*(shift-delay))
        end
    end
    for iw = nw+1:nf
        D[iw, :] = conj(D[nf-iw+2, :])
    end
    d = ifft(D, 1)
    d = real(d[1:nt, :])
    extent = Extent(Int32(nt), Int32(nx), Int32(1), Int32(1), Int32(1),
                    Float32(ot), Float32(ox), Float32(0), Float32(0),
                    Float32(0), Float32(dt), Float32(dx), Float32(1),
                    Float32(1), Float32(1), "Time", "Offset", "ix2", "ix3",
                    "ix4", "s", "m", "index", "index", "index",
                    "Hyperbolic events")
    return d, extent

end
