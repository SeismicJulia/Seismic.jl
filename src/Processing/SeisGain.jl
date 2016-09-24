"""
    SeisGain(d, t; <keyword arguments>)

Gain a group of traces.

# Arguments
* `d::Array{Real,2}`: two dimensional data.

# Keyword arguments
* `dt::Real=0.002`: sampling interval in secs.
* `kind::String="time"`: if kind="time", gain = t.^a . * exp(-bt);
                         if kind="agc", automatic gain control is applied.
* `param::Vector{Real}=[2.0,0.0]`: if kind="time", param = [a,b];
                                   if kind="agc", param = [agc_gate]
* `norm::Int=0`: `norm=0` no normalization; `norm=1` normalize each trace by
                  amplitude; `norm=2` normalize each trace by rms value/

# Output
* d::Array{Real, 2}`: gained two dimensional data.

# Example
```julia
julia> d, extent = SeisHypEvents(); 
       dout = SeisGain(d, kind="agc", param=[0.05]); 
       SeisPlot([d dout], extent);
```

Credits: Juan I. Sabbione, Aaron Staton, Mauricio D. Sacchi, 2016

"""
function SeisGain{Td<:Real,Tp<:Real
                  }(d::Array{Td,2}; dt::Real=0.004, kind::String="time",
                    param::Vector{Tp}=[2.0,0.0], norm::Int=0)
    
    nt = size(d,1)
    nx = size(d[:,:],2)
    dout = zeros(d)

    if kind == "time"   # Geometrical spreading-like gain
        
        a = param[1]
        b = param[2]
        t = collect(0:1:nt-1)*dt
        tgain = (t.^a).*exp(b.*t)
        for k = 1:nx
            dout[:,k] = d[:,k].*tgain
        end
        
    end

    if kind=="agc"   # AGC 
        
        L = floor(Int,round(param[1]/(2dt)))
        h = triang(2*L+1)
        
        for k = 1:nx
            aux =  d[:,k]
            e = aux.^2
            rms = sqrt(abs(conv(e,h)[L+1:nt+L]))
            epsi = 1.e-10*maximum(rms)
            op = rms./(rms.^2+epsi)
            dout[:,k] = d[:,k].*op
        end
    end
    
    if norm==1     # Normalize by amplitude 
        
        for k = 1:nx
            aux =  d[:,k]
            amax = maximum(abs(aux))
            dout[:,k] = dout[:,k]/amax
        end
        
    end
    
    if norm==2;    # Normalize by rms 
        
        for k = 1:nx
            aux = d[:,k];
            amax = sqrt(sum(aux.^2)/nt);
            dout[:,k] = dout[:,k]/amax;
        end
        
    end

    return dout

end

function triang(n::Integer)
    [1 - abs((k - (n-1)/2))/(n/2) for k=0:(n-1)]
end
