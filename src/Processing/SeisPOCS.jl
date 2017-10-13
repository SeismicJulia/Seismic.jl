"""
**SeisPOCS**

*Projection Onto Convex Sets interpolation of seismic records.*

**IN**

* d_in: input data that can have up to 5 dimensions
* p=1., exponent for thresholding (1 is equivalent to soft thres. high number is equivalent to hard thresholding)
* alpha=1 add-back ratio for imputation step. Use 1 for noise free data, and < 1 for denoising of original traces.
* dt=0.001 sampling rate along the time axis (in seconds)
* fmax=99999. maximum temporal frequency to process.
* padt=2 padding to use for the time axis
* padx=1 padding to use for the spatial axes
* Niter=100 number of iterations

**OUT**

* d_out: interpolated data
"""
function SeisPOCS(in;p=1.,dt=0.001,fmax=99999.,padt=2,padx=1,Niter=100,alpha=1)

    perci = 0.999999;
    percf = 0.0;
    nt = size(in,1)
    nx1 = size(in,2)
    nx2 = size(in,3)
    nx3 = size(in,4)
    nx4 = size(in,5)
    d = zeros(Float32,nt,nx1,nx2,nx3,nx4)
    d[1:nt,1:nx1,1:nx2,1:nx3,1:nx4] = in
    nf = padt*nextpow2(nt)
    dw = 2.*pi/nf/dt
    nw = round(Int,nf/2) + 1
    fmax = fmax < 0.5/dt ? fmax : 0.5/dt
    if(fmax*dt*nf < nw)
	iw_max = round(Int,floor(fmax*dt*nf))
    else
	iw_max = round(Int,floor(0.5/dt))
    end
    nx1 > 1 ? nk1 = padx*nextpow2(nx1) : nk1 = 1
    nx2 > 1 ? nk2 = padx*nextpow2(nx2) : nk2 = 1
    nx3 > 1 ? nk3 = padx*nextpow2(nx3) : nk3 = 1
    nx4 > 1 ? nk4 = padx*nextpow2(nx4) : nk4 = 1
    nk = nk1*nk2*nk3*nk4
    # generate sampling operator from the padded data
    T = CalculateSampling(d)
    T = T[1:1,:,:,:,:]
    if (sum(T[:])/length(T[:]) < 0.05)
	println(sum(T[:])/length(T[:]))
	return in
    else
	T = Pad5D(T,1,nk1,nk2,nk3,nk4)
	T = T[1,:,:,:,:]
	d = Pad5D(d,nf,nk1,nk2,nk3,nk4)
	D = fft(d,1)
	for iw=1:iw_max
	    x = D[iw,:,:,:,:]
	    y = copy(x)
	    for iter = 1 : Niter
		Y = fft(y)
		amp = sort(vec(abs.(Y[:])))
		perc = perci + (iter-1)*((percf-perci)/Niter);
		cutoff = amp[round(Int,floor(perc*nk)+1)];
		for j = 1 : nk1*nk2*nk3*nk4
		    if (abs(Y[j]) < cutoff)
			Y[j] = 0.
		    else
			Y[j] = Y[j]*(1 - (cutoff/(abs(Y[j]) + 0.0000001))^p)
		    end
		end
		y = ifft(Y)
		y = alpha*x + (1-alpha*T).*y
	    end
	    D[iw,:,:,:,:] = y
	end
	# symmetries
	for iw=nw+1:nf
	    D[iw,:,:,:,:] = conj(D[nf-iw+2,:,:,:,:])
	end
	d = ifft(D,1)
	d = real(d[1:nt,1:nx1,1:nx2,1:nx3,1:nx4])
	return d
    end
end
