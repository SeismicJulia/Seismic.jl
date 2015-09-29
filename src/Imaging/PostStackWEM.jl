function PostStackWEM(m,d,param)
	# 
	# PostStackWEM: Post Stack Wave Equation Migration and Demigration 
	# of 3D isotropic data.

	adj = get(param,"adj",true) # flag for adjoint (migration), or forward (demigration)
	nt = get(param,"nt",1) # trace length
	ot = get(param,"ot",0) # trace origin
	dt = get(param,"dt",0.001) # trace sample rate
	damping = get(param,"damping",1000.)  # damping for deconvolution imaging condition
	vel = get(param,"vel","vel") # seis file containing the velocity (should have same x and z dimensions as the desired image)
	wav = get(param,"wav","wav") # seis file containing the source wavelet (in time domain) 
	gz = get(param,"gz",0.) # receiver depth
	pade_flag = get(param,"pade_flag",false) # flag for Pade Fourier correction
	fmin = get(param,"fmin",0)  # min frequency to process (Hz)
	fmax = get(param,"fmax",floor(0.5/dt))  # min frequency to process (Hz)
	padt = get(param,"padt",2)  # pad factor for the time axis
	padx = get(param,"padx",2) # pad factor for the spatial axes
	omp = get(param,"omp",1) # number of shared memory threads to use for frequency slice processing 
	verbose = get(param,"verbose",false) # flag for error / debugging messages

	v,h = SeisRead(vel)
	min_imx = h[1].imx
	max_imx = h[end].imx
	nmx = max_imx - min_imx + 1
	min_imy = h[1].imy
	max_imy = h[end].imy
	nmy = max_imy - min_imy + 1
	dmx = h[nmy+1].mx - h[nmy].mx > 0 ? h[nmy+1].mx - h[nmy].mx : 1.
	dmy = h[2].my - h[1].my > 0 ? h[2].my - h[1].my : 1.
	nz = int(h[1].n1)
	dz = h[1].d1
	oz = h[1].o1
	param_zowem = {"adj"=>adj,
	               "omx"=>min_imx,"dmx"=>dmx,"nmx"=>nmx,
	               "omy"=>min_imy,"dmy"=>dmy,"nmy"=>nmy,
	               "oz"=>oz,"dz"=>dz,"nz"=>nz,
	               "ot"=>ot,"dt"=>dt,"nt"=>nt,	               
	               "padt"=>padt,"padx"=>padx,
	               "v"=>v,
	               "fmin"=>fmin,"fmax"=>fmax}
	if (adj)
		d1,h = SeisRead(d)
		d1 = reshape(d1,nt,nmx,nmy)
		m1 = zowem(d1,param_zowem)
		for ix = 1 : nmx*nmy
			h[ix].n1 = nz
			h[ix].o1 = oz
			h[ix].d1 = dz
		end
		SeisWrite(m,m1,h)
	else
		m1,h = SeisRead(m)
		m1 = reshape(m1,nz,nmx,nmy)
		d1 = zowem(m1,param_zowem)
		for ix = 1 : nmx*nmy
			h[ix].n1 = nt
			h[ix].o1 = ot
			h[ix].d1 = dt
		end
		SeisWrite(d,d1,h)
	end
	
	return
end

function zowem(in,param::Dict{Any,Any})

	adj = get(param,"adj",true)
	omx = get(param,"omx",0)
	dmx = get(param,"dmx",1)
	nmx = get(param,"nmx",1)
	omy = get(param,"omy",0)
	dmy = get(param,"dmy",1)
	nmy = get(param,"nmy",1)
	oz = get(param,"oz",0)
	dz = get(param,"dz",1)
	nz = get(param,"nz",1)
	ot = get(param,"ot",0)
	dt = get(param,"dt",1)
	nt = get(param,"nt",1)
	padt = get(param,"padt",2)
	padx = get(param,"padx",2)
	v = get(param,"v",[])
	fmin = get(param,"fmin",0)
	fmax = get(param,"fmax",floor(0.5/dt))
	
	nf = padt*nt
	nkx = nmx > 1 ? padx*nmx : 1
	nky = nmy > 1 ? pady*nmy : 1
	dkx = 2*pi/nkx/dmx;
	dky = 2*pi/nky/dmy;
	dw = 2.*pi/nf/dt
	nw = int(nf/2) + 1
	if (adj)
		d = pad3d(in,nf,nkx,nky)
		m = zeros(Float32,nz,nkx,nky)
	else
		m = pad3d(in,nz,nkx,nky)
		d = zeros(Float32,nf,nkx,nky)
	end	

	# decompose slowness into layer average, and layer perturbation
	po = zeros(Float32,nz) 
	pd = zeros(Float32,nz,nmx*nmy)
	for iz = 1 : nz
		po[iz] = 1./mean(v[iz,:])
		pd[iz,:] = 1./v[iz,:] - po[iz]
	end
	pd = reshape(pd,nz,nmx,nmy)
	pd = pad3d(pd,nz,nkx,nky)
	param["po"] = po
	param["pd"] = pd
	
	if(1 <= int(fmin*dt*nf) < nw) 
		iw_min = int(floor(fmin*dt*nf))
	else 
		iw_min = 1
	end
	if(fmax*dt*nf < nw) 
		iw_max = int(floor(fmax*dt*nf))
	else 
		iw_max = int(floor(0.5/dt))
	end
	# set up kx and ky vectors
	kx = zeros(Float32,nkx,nky)	
	ky = zeros(Float32,nkx,nky)
	for ikx = 1 : nkx
		kx[ikx,:] = ikx <= int(nkx/2) ? dkx*(ikx-1) : -(dkx*(nkx-1) - dkx*(ikx-1))
	end
	for iky = 1 : nky
		ky[:,iky] = iky <= int(nky/2) ? dky*(iky-1) : -(dky*(nky-1) - dky*(iky-1))
	end
	param["kx"] = kx
	param["ky"] = ky
	if (adj)
		D = fft(d,1)
		for iw=iw_min:iw_max
			param["w"] = iw*dw
			x = squeeze(D[iw,:,:],1)
			m,x = extrap1f(m,x,param)
			D[iw,:,:] = x
		end
		return real(m[1:nz,1:nmx,1:nmy])
	else
		D = zeros(Complex{Float32},nf,nkx,nky)
		for iw=iw_min:iw_max
			param["w"] = iw*dw
			x = squeeze(D[iw,:,:],1)
			m,x = extrap1f(m,x,param)
			D[iw,:,:] = x
		end
		# symmetries
		for iw=nw+1:nf
			D[iw,:,:] = conj(D[nf-iw+2,:,:])
		end 
		d = ifft(D,1)
		return real(d[1:nt,1:nmx,1:nmy])
	end	
	
end

function extrap1f(m,d,param)

	adj = get(param,"adj",true)
	w = get(param,"w",1)
	dz = get(param,"dz",1)
	po = get(param,"po",[])
	pd = get(param,"pd",[])
	nz = get(param,"nz",1)
	kx = get(param,"kx",[])
	ky = get(param,"ky",[])

	if (adj)
		for iz = 1 : 1 : nz
			d = fft_op(d,{"adj"=>true})
			s = (w.^2)*(po[iz]*po[iz]) - kx.^2 - ky.^2
			s[(s.<.0)] = 0.
			d = d.*exp(1im*sqrt(s)*dz)
			d = fft_op(d,{"adj"=>false})
			#d = d.*exp(1im*w*pd[iz]*dz)
			m[iz,:,:] = squeeze(m[iz,:,:],1) + 2*real(d)
		end
	else
		for iz = nz : -1 : 1
			d += squeeze(m[iz,:,:],1)
			#d = d.*exp(-1im*w*pd[iz]*dz)
			d = fft_op(d,{"adj"=>true})
			s = (w.^2)*(po[iz]^2) - kx.^2 - ky.^2
			s[(s.<.0)] = 0.
			d = d.*exp(-1im*sqrt(s)*dz)	
			d = fft_op(d,{"adj"=>false})
		end
	end
	return m,d
	
end

function pad3d(a,N1,N2,N3)
	n1,n2,n3 = size(a)
	b = zeros(Float32,N1,N2,N3)
	b[1:n1,1:n2,1:n3] = a
	return b
end
