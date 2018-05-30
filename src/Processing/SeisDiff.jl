"""
```
SeisDiff(d ; <keyword arguments>)
```

Apply differentiation to seismic traces in freq. It can also be 
used to apply a phase rotation 

# Arguments
* `d`: Input 2D data array in tx domain. First dimension is time.

# Keyword arguments
* `delay=0.1`: desired time delay in seconds.
* `pow=-2`: order of derivative
* `rot=0`: constant phase shift or rotation 


# Example
```
julia> d,extent = SeisLinearEvents(d); SeisPlot(d);
julia> d2 = SeisDiff(d;dt=0.004,pow=1); SeisPlot(d);
```
"""


function SeisDiff(d;dt=0.001,pow=-2,rot=0)


	nt = size(d,1)
	D = fft(d,1)
	dw = 2.*pi/nt/dt
	nw = Int(nt/2) + 1
	eps = pow < 0 ? 10.0 : 0.0
	for iw=1:nw
		D[iw,:] *= exp(rot*1im)*((iw*dw) + eps)^pow
	end
	# symmetries
	for iw=nw+1:nt
		D[iw,:] = conj(D[nt-iw+2,:])
	end
	d = real(ifft(D,1))
	return d

end

function SeisDiff(in,h::Array{Header,1};pow=-2,rot=0)

	out = differentiate_traces(in;dt=h[1].d1,pow=pow,rot=rot)
	return out,h

end

function SeisDiff(in::AbstractString,out::AbstractString;pow=-2,rot=0)

	@compat parameters = Dict(:pow=>pow,:rot=>rot)
	SeisProcess(in,out,[SeisDiff],[parameters])

end
