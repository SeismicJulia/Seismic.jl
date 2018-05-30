"""
```
SeisEnvelope(d)
```

Calculate the envelope attribute of  a group of input traces

# Arguments
* `d`: Input data. A 2D array where the first dimension is time.

# Examples
```
julia> dtsec = 0.002; w = Ricker(dt=dtsec); 
julia> e = SeisEnvelope(w); t=dtsec*collect(0:1:length(w)-1);plot(t,w,t,e,"-r");xlabel("Time [s]")

julia> d,extent = SeisLinearEvents(); SeisPlot(d,extent,style="wiggles")
julia> out = SeisEnvelope(d); SeisPlot(out,extent,style="wiggles")
```
"""

function SeisEnvelope(d)

	D = fft(d,1)
	D[1:Int(floor(size(d,1)/2)),:] = 0.0
	return 2*abs.(ifft(D,1))

end

function SeisEnvelope(A,h::Array{Header,1})
	for itrace = 1 : size(A[:,:],2)
		A[:,itrace] = SeisEnvelope(A[:,itrace])
	end
	return A,h
end

function SeisEnvelope(in::AbstractString,out::AbstractString)
	@compat parameters = Dict()
	SeisProcess(in,out,[SeisEnvelope],[parameters];group="some",ntrace=100000)
end
