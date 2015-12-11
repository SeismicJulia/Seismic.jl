function SeisBandPass(d,h::Array{Header,1},param=Dict())

	param["ot"] = h[1].o1
	param["dt"] = h[1].d1
	d = SeisBandPass(d,param)

	return d,h;
end

function SeisBandPass(d,param::Dict{Any,Any})

	nt = size(d,1)
	nx = size(d[:,:],2)
	ot = get(param,"ot",0)
	dt = get(param,"dt",0.001)
	nf = iseven(nt) ? nt : nt + 1
	df = 1/nf/dt
	nw = int(nf/2) + 1
	fmin = get(param,"fmin",0)
	fmax = get(param,"fmax",0.5/dt)

	fa = fmin < 5 ? 0. : fmin - 5.
	fb = fmin
	fc = fmax	
	fd = fmax > 0.5/dt - 5. ? 0.5/dt : fmax + 5.

	if(fd*dt*nf < nw) 
		iw_max = int(floor(fd*dt*nf))
	else 
		iw_max = int(floor(0.5/dt))
	end

	d = pad_first_axis(d,nf)
	m = fft(d,1)/sqrt(size(d,1))

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
	#println("size(m)= ",size(m))
	#println("iw_max= ",iw_max)
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
