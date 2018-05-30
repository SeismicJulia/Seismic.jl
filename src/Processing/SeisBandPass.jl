"""
```
SeisBandPass(d ; <keyword arguments>)
```

Apply a bandpass filter to seismic data. Input and output data is in tx domain. Filter is applied in fx domain.

# Arguments
* `d`: Input 2D data array in tx domain. Time is first dimension. 

# Keyword arguments
* `dt=0.001`: time sampling interval.
* `fa=0,fb=0,fc=60,fd=80`: corner frequencies in Hz. 


# Output
* `d`: Filtered 2d data array in tx domain.

# Examples
```
julia> d,e = SeisLinearEvents(); SeisPlot(d,plot_type="Amplitude",dy=0.004,fmax=125);
julia> d_filter = SeisBandPass(d;dt=0.004,fa=2,fb=8,fc=12,fd=20); SeisPlot(d_filter,plot_type="Amplitude",dy=0.004,fmax=125);
julia> SeisPlot([d d_filter],e,title="Data and Filtered data")



```
"""

function SeisBandPass(d;dt=0.001,fa=0,fb=0,fc=60,fd=80)

	nt = size(d,1)
	nx = size(d[:,:],2)
	nf = iseven(nt) ? nt : nt + 1
	df = 1/nf/dt
	nw = round(Int,nf/2) + 1

	if(fd*dt*nf < nw)
		iw_max = round(Int,floor(fd*dt*nf))
	else
		iw_max = round(Int,floor(0.5/dt))
	end

	d = pad_first_axis(d,nf)
	m = fft(d,1)/sqrt(size(d,1))
	if fa > 0.
		m[1,:] *= 0.
	end
	for iw=2:iw_max
		f = df*(iw-1)
		if (f<fa)
			m[iw,:] *= 0.
		elseif (f >= fa && f < fb)
			m[iw,:] *= (f-fa)/(fb-fa)
		elseif (f >= fb && f <= fc)
			m[iw,:] *= 1.
		elseif (f > fc && f <= fd)
			m[iw,:] *= 1. - (f-fc)/(fd-fc)
		else
			m[iw,:] *= 0.
		end
	end
	m[iw_max:end,:] = 0.

	# symmetries
	for iw=nw+1:nf
		m[iw,:] = conj(m[nf-iw+2,:])
	end
	d = real(bfft(m,1)/sqrt(size(m,1)))

	return d[1:nt,1:nx];
end

function pad_first_axis(a,N1)
	n1 = size(a,1)
	nx = size(a[:,:],2)
	b = zeros(N1,nx)
	for ix = 1 : nx
		b[1:n1,ix] = a[:,ix]
	end
	return b
end

function SeisBandPass(d,h::Array{Header,1};dt=0.001,fa=0,fb=0,fc=60,fd=80)

	out = SeisBandPass(d;dt=h[1].d1,fa=fa,fb=fb,fc=fc,fd=fd)
	return out,h

end

function SeisBandPass(d,h::Array{Header,1},ext;dt=0.001,fa=0,fb=0,fc=60,fd=80)

	out = SeisBandPass(d;dt=h[1].d1,fa=fa,fb=fb,fc=fc,fd=fd)
	return out,h,ext

end


function SeisBandPass(in::AbstractString,out::AbstractString;fa=0,fb=0,fc=60,fd=80)

	@compat parameters = Dict(:fa=>fa,:fb=>fb,:fc=>fc,:fd=>fd)
	SeisProcess(in,out,[SeisBandPass],[parameters])

end
