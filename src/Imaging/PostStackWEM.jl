function PostStackWEM(m,d,adj;vel="NULL",damping=1000,nt=1,ot=0,dt=0.001,wav="NULL",gz=0,fmin=0,fmax=99999,padt=2,padx=2,verbose=true)
    # 
    # PostStackWEM: Post Stack Wave Equation Migration and Demigration 
    # of 3D isotropic data.
    
    #adj: flag for adjoint (migration), or forward (demigration)
    #nt: trace length
    #ot:trace origin
    #dt : trace sample rate
    #damping: damping for deconvolution imaging condition
    #vel :seis file containing the velocity (should have same x and z dimensions as the desired image)
    #wav :seis file containing the source wavelet (in time domain) 
    #gz : receiver depth
    #fmin : min frequency to process (Hz)
    #fmax : max frequency to process (Hz)
    #padt : pad factor for the time axis
    #padx : pad factor for the spatial axes 
    #verbose : flag for error / debugging messages
    
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
	d1 = zowem(m1,adj,omx=min_imx,dmx=dmx,nmx=nmx,omy=min_imy,dmy=dmy,nmy=nmy,oz=oz,dz=dz,nz=nz,ot=ot,dt=dt,nt=nt,padt=padt,padx=padx,v=v,fmin=fmin,fmax=fmax)
	for ix = 1 : nmx*nmy
	    h[ix].n1 = nt
	    h[ix].o1 = ot
	    h[ix].d1 = dt
	end
	SeisWrite(d,d1,h)
    end
    
    return
end

function zowem(in,adj;omx=0,dmx=1,nmx=1,omy=0,dmy=1,nmy=1,oz=0,dz=1,nz=1,ot=0,dt=1,nt=1,padt=2,padx=2,v=[],fmin=0.,fmax=999999.)
    
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
	    x = D[iw,:,:]
	    m,x = extrap1f(m,x,param)
	    D[iw,:,:] = x
	end
	return real(m[1:nz,1:nmx,1:nmy])
	else
	    D = zeros(Complex{Float32},nf,nkx,nky)
	for iw=iw_min:iw_max
	    param["w"] = iw*dw
	    x = D[iw,:,:]
	    m,x = extrap1f(m,x,adj,w=w,dz=dz,po=po,pd=pd,nz=nz,kx=kx,ky=ky)
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

function extrap1f(m,d,adj,w=1,dz=1,po=[],pd=[],nz=1,kx=[],ky=[])

    if (adj)
	for iz = 1 : 1 : nz
	    d = FFTOp(d,true)
	    s = (w.^2)*(po[iz]*po[iz]) - kx.^2 - ky.^2
	    s[(s.<.0)] = 0.
	    d = d.*exp(1im*sqrt(s)*dz)
	    d = fft_op(d,false)
	    #d = d.*exp(1im*w*pd[iz]*dz)
	    m[iz,:,:] = m[iz,:,:] + 2*real(d)
	end
    else
	for iz = nz : -1 : 1
	    d += m[iz,:,:]
	    #d = d.*exp(-1im*w*pd[iz]*dz)
	    d = fft_op(d,true)
	    s = (w.^2)*(po[iz]^2) - kx.^2 - ky.^2
	    s[(s.<.0)] = 0.
	    d = d.*exp(-1im*sqrt(s)*dz)	
	    d = fft_op(d,false)
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
