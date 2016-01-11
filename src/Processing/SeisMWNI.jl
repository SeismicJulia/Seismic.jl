"""
**SeisMWNI**

*Minimum Weighted Norm Interpolation of seismic records.*

**IN**   

* d_in: input data that can have up to 5 dimensions
* dt=0.001 sampling rate along the time axis (in seconds)
* fmax=99999. maximum temporal frequency to process.
* padt=2 padding to use for the time axis
* padx=1 padding to use for the spatial axes
* Niter_internal=10 number of internal iterations for Conjugate Gradients
* Niter_external=3 number of external iterations for iterative reweighting

**OUT**  

* d_out: interpolated data
"""
function SeisMWNI(in,dt=0.001,fmax=99999.,padt=2,padx=1,Niter_internal=10,Niter_external=3)

    nt = size(in,1)
	nx1 = size(in,2)
	nx2 = size(in,3)
	nx3 = size(in,4)
	nx4 = size(in,5)
	d = zeros(Float32,nt,nx1,nx2,nx3,nx4)
	d[1:nt,1:nx1,1:nx2,1:nx3,1:nx4] = in
	nf = padt*nextpow2(nt)
	dt = get(param,"dt",0.001)
	dw = 2.*pi/nf/dt
	nw = int(nf/2) + 1
	fmax = get(param,"fmax",int(floor(0.5/dt)))
	if(fmax*dt*nf < nw) 
		iw_max = int(floor(fmax*dt*nf))
	else 
		iw_max = int(floor(0.5/dt))
	end
	nx1 > 1 ? nk1 = padx*nextpow2(nx1) : nk1 = 1
	nx2 > 1 ? nk2 = padx*nextpow2(nx2) : nk2 = 1
	nx3 > 1 ? nk3 = padx*nextpow2(nx3) : nk3 = 1
	nx4 > 1 ? nk4 = padx*nextpow2(nx4) : nk4 = 1
	nk = nk1*nk2*nk3*nk4
	# generate sampling operator from the padded data
	T = CalculateSampling(d)
	T = T[1,:,:,:,:]
	param_dataweights = Dict(:wd=>T) 
    param_fft_op =  Dict(:wd=>T) 
    param_modelweights = Dict(:wm=>wm)

	if (sum(T[:])/length(T[:]) < 0.05)
		println(sum(T[:])/length(T[:]))
		return in
	else
		T = Pad5D(T,1,nk1,nk2,nk3,nk4)
		T = squeeze(T[1,:,:,:,:],1)
		param["wd"] = T	
		d = Pad5D(d,nf,nk1,nk2,nk3,nk4)
		D = fft(d,1)
		for iw=1:iw_max
			x = squeeze(D[iw,:,:,:,:],1)
			y = copy(x)
			Y = fft_op(y,adj=true)
			Y,misfit = IRLS(Y.*0,x,Niter_external,Niter_internal,[ApplyDataWeights fft_op ApplyModelWeights],[param_dataweights param_fft_op param_modelweights])
			y = fft_op(Y,adj=false)
			D[iw,:,:,:,:] = y
		end
		# symmetries
		for iw=nw+1:nf
			D[iw,:,:,:,:] = conj(D[nf-iw+2,:,:,:,:])
		end 
		d = ifft(D,1)
		d = real(d[1:nt,1:nx1,1:nx2,1:nx3,1:nx4])
		return reshape(d,size(in))
	end
end
