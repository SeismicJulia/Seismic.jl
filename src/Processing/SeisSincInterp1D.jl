###
"""
```
SeisSincInterp1D(d,order)
```

 Resample seismic traces via 1D sinc interpolation. The time series are sampled every dt secs, 
 the output corresponds to series with time interval dt/`order`.

 # Arguments
 * `d::Array{Real,N}`: N-dimensional data, first dimension is time
 * `order::Integer`: order of interpolation 2,4

 # Examples

```
# 2-time upsampling of a wavelet 
julia> order = 2; dt=0.004; w = Ricker(dt=dt, f0=20); t = collect(0:1:length(w)-1)*dt;
julia> wout = SeisSincInterp1D(w,order); tout = collect(0:1:length(wout)-1)*dt/order;
julia> plot(t,w); plot(tout,wout,"*")  


# 4-time upsampling of a gather 
julia> d,e = SeisLinearEvents(); di = SeisSincInterp1D(d,4);
julia> figure(1); clf(); SeisPlot(d,style="wiggles")
julia> figure(2); clf(); SeisPlot(di,style="wiggles")  


```
MAY 2018, MDS

"""


 function SeisSincInterp1D(d::Array{Tv},order::Ti) where {Tv<:Number, Ti<:Int64} 


        Ndims = ndims(d)
        Ndims == 1? nx = 1 :  nx = prod(size(d)[2:Ndims])
        dims = size(d) 
        N = dims[1]
        Npad = order*N 
        nf = 2*Npad
        nw = convert(Int,floor(nf/2)) + 1
        di = zeros(eltype(d),nf,dims[2:end]...)

        for k = 1:nx

        dd = zeros(eltype(d),nf)
        dd[1:order:Npad] = d[1:1:end,k]
        D = fft(dd)
        nyq = convert(Int,floor(nw/order)) + 1
        for iw = nyq : nw
                D[iw] *= 0
        end
        # symmetries
        for iw=nw+1:nf
                D[iw] = conj(D[nf-iw+2])
        end
        di[:,k] = real(ifft(D,1))
        end

        return  reshape(di[1:Npad,:]*order,Npad,dims[2:end]...)
        

end
