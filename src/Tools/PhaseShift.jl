function PhaseShift(x,angle)

	nx = length(x)
	nf = 2*nx
	nw = Int(nf/2) + 1
	X = fft([x;zeros(typeof(x[1]),nf-nx)])
	for iw=1:nw
		X[iw] *= exp(-1im*angle)
	end
	# symmetries
	for iw=nw+1:nf
		X[iw] = conj(X[nf-iw+2])
	end
	x = real(ifft(X,1))
	return x[1:nx]

end
