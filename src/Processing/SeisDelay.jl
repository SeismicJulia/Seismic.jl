function SeisDelay(d; delay=0.1, dt=0.001)
	
    nt = size(d,1)
    nf = nt
    D = fft(d,1);
    dw = 2*pi/nt/dt
    nw = round(Int,nt/2) + 1
    for iw = 1 : nw
	w = (iw-1)*dw
	D[iw,:] = D[iw,:]*exp(-1im*w*delay)
    end
    # symmetries
    for iw=nw+1:nf
	D[iw,:] = conj(D[nf-iw+2,:])
    end 
    d2 = ifft(D,1)
    d2 = real(d2[1:nt,:])
    return d2
    
end
