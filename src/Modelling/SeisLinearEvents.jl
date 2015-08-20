function Ricker(f0,dt)
	nw=int(2.2/f0/dt)
	nw=2*int(nw/2)+1
	nc=int(nw/2)
	w = zeros(nw,1)
	k=[1:1:nw]'
	alpha = (nc-k+1)*f0*dt*pi
	beta=alpha.^2
	w = (1.-beta.*2).*exp(-beta)
	return w
end

function SeisLinearEvents(param=Dict())	
	ot = get(param,"ot",0.)
	dt = get(param,"dt",0.004)
	nt = get(param,"nt",500)
	ox1 = get(param,"ox1",0.)
	dx1 = get(param,"dx1",10)
	nx1 = get(param,"nx1",20)
	ox2 = get(param,"ox2",0.)
	dx2 = get(param,"dx2",10)
	nx2 = get(param,"nx2",1)
	ox3 = get(param,"ox3",0.)
	dx3 = get(param,"dx3",10)
	nx3 = get(param,"nx3",1)
	ox4 = get(param,"ox4",0.)
	dx4 = get(param,"dx4",10)
	nx4 = get(param,"nx4",1)
	tau1 = get(param,"tau1",0.3)
	tau2 = get(param,"tau2",0.)
	tau3 = get(param,"tau3",0.)
	tau4 = get(param,"tau4",0.)
	v1 = get(param,"v1",1500.)
	v2 = get(param,"v2",1500.)
	v3 = get(param,"v3",1500.)
	v4 = get(param,"v4",1500.)
	amp = get(param,"amp",1.)
	f0 = get(param,"f0",20.)

	d = zeros(nt,nx1,nx2,nx3,nx4)	
	nf = 4*nextpow2(nt)
	dw = 2.*pi/nf/dt
	nw = int(nf/2) + 1
	t = [0:1:nt-1]'*dt
	x1 = Array(Float32,1,nx1,nx2,nx3,nx4)
	x2 = Array(Float32,1,nx1,nx2,nx3,nx4)
	x3 = Array(Float32,1,nx1,nx2,nx3,nx4)
	x4 = Array(Float32,1,nx1,nx2,nx3,nx4)
    for ix1=1:nx1 
    	for ix2=1:nx2 
    		for ix3=1:nx3 
    			for ix4=1:nx4  
    				x1[1,ix1,ix2,ix3,ix4] = ix1*dx1 + ox1
    				x2[1,ix1,ix2,ix3,ix4] = ix2*dx2 + ox2
    				x3[1,ix1,ix2,ix3,ix4] = ix3*dx3 + ox3
    				x4[1,ix1,ix2,ix3,ix4] = ix4*dx4 + ox4
    			end
    		end
    	end
    end

	D = zeros(Complex{Float64},nf,nx1,nx2,nx3,nx4)
	for ievent=1:length(tau1)
		wav = Ricker(f0[ievent],dt)
    	nwav = length(wav)
    	wav = cat(2,wav,zeros(1,nf-length(wav)))
    	Wav = fft(wav',1)
		delay = dt*(int(nwav/2)+1)
		for iw=1:nw
    		w = (iw-1)*dw
    		shift = tau1[ievent] + tau2[ievent] + tau3[ievent] + tau4[ievent] + x1/v1[ievent] + x2/v2[ievent] + x3/v3[ievent] + x4/v4[ievent] - delay
        	D[iw,:,:,:,:] += amp[ievent]*Wav[iw]*exp(-1im*w*shift)
    	end
    end
    # symmetries
    for iw=nw+1:nf
		D[iw,:,:,:,:] = conj(D[nf-iw+2,:,:,:,:])
    end 
	d = ifft(D,1)
	d = real(d[1:nt,:,:,:,:])
	return d
end
