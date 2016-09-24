function differentiate_traces(d;dt=0.001,pow=-2,rot=0)

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

function SeisDiff(in::String,out::String;pow=-2,rot=0)

	@compat parameters = Dict(:pow=>pow,:rot=>rot)
	SeisProcess(in,out,[SeisDiff],[parameters])

end
