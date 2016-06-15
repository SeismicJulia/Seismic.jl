"""
    SeisAddNoise(d, snr; <keyword arguments>)

Add noise at a given signal-to-noise ratio level `snr` to a N-dimensional input
data `d`. Noise can be band limited using kewyord `L`.

# Arguments
* `d::Array{Real, N}`: N-dimensional data. `N` must be <= 5.
* `snr::Real`: signal-to-noise ratio.

# Keyword arguments
* `db::Bool=false`: `db=false` if `snr` is given by amplitude, `db=false` if
snr is given in dB.
* `pdf::ASCIIString="gaussian"`: random noise probability distribution:
`"gaussian"` or `"uniform"`.
* `L::Int=1`: averaging operator length to band-limit the random noise.

# Examples
```julia
julia> w = Ricker(); wn = SeisAddNoise(w, 2); plot(w); plot(wn); SNR(w, wn)

julia> d, extent = SeisHypEvents(); dn = SeisAddNoise(d, 1.0, db=true, L=9);
SeisPlot([d dn], extent); SNR(d, dn, db=true)
```
"""

function SeisAddNoise{T<:Real, N}(d::Array{T, N}, snr::Real; db::Bool=false, 
                                  pdf::ASCIIString="gaussian", L::Int=1)

    noise = GenNoise(size(d), pdf, L=L)
    if db==false
        noise = noise/norm(noise) * norm(d)/snr
    elseif db==true
        noise = noise/norm(noise) * norm(d)/10.0^(0.05*snr)
    end
    noisy = d + noise
    assert(abs(SNR(d, noisy, db=db) - snr) < 1e10*eps(AbstractFloat(snr)))
    return noisy

end

# Generate the noise 
function GenNoise(dims::Tuple, pdf::ASCIIString; L::Int=1)

    N = length(dims)
    n1 = dims[1]
    hamm = Hamming(L)
    L2 = floor(Int,L/2)
    if pdf=="gaussian"
        noise = randn(dims)
    elseif pdf=="uniform"
        noise = -1 + 2*rand(dims)
    else
        error("pdf must be gussian or uniform")
    end
    if N==1
        noise = conv(noise, hamm)[L2+1:n1+L2]
    elseif N==2
        for i2 = 1:dims[2]
            noise[:,i2] = conv(noise[:,i2], hamm)[L2+1:n1+L2]
        end
    elseif N==3
        for i2 = 1:dims[2], i3 = 1:dims[3]
            noise[:,i2,i3] = conv(noise[:,i2,i3], hamm)[L2+1:n1+L2]
        end
    elseif N==4
        for i2 = 1:dims[2], i3 = 1:dims[3], i4 = 1:dims[4]
            noise[:,i2,i3,i4] = conv(noise[:,i2,i3,i4], hamm)[L2+1:n1+L2]
        end
    elseif N==5
        for i2 = 1:dims[2], i3 = 1:dims[3], i4 = 1:dims[4], i5 = 1:dims[5]
            noise[:,i2,i3,i4,i5] = conv(noise[:,i2,i3,i4,i5], hamm)[L2+1:n1+L2]
        end
    else
        error("Maximum dimension of array is N=5")
    end
    return noise

end
