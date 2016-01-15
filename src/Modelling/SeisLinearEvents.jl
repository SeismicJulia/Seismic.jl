"""
**SeisLinearEvents**

*Generate five dimensional data consisting of linear events*

**IN**   

* ot=0
* dt=0.004
* nt=500
* ox1=0
* dx1=10
* nx1=100
* ox2=0
* dx2=10
* nx2=1
* ox3=0
* dx3=10
* nx3=1
* ox4=0
* dx4=10
* nx4=1
* tau1=[1.0,1.6]
* tau2=0
* tau3=0
* tau4=0
* v1=[2500,-800]
* v2=1500
* v3=1500
* v4=1500
* amp=[1,-1]
* f0=20
* ricker=true
* exponent=1
* sinusoidal=false

**OUT**  

* d: modelled data with dimensions d[1:nt,1:nx1,1:nx2,1:nx3,1:nx4]

*Credits: Aaron Stanton, 2015*

"""
function SeisLinearEvents(;ot=0,dt=0.004,nt=500,ox1=0,dx1=10,nx1=100,ox2=0,dx2=10,nx2=1,ox3=0,dx3=10,nx3=1,ox4=0,dx4=10,nx4=1,tau1=[1.0,1.6],tau2=0,tau3=0,tau4=0,v1=[2500,-800],v2=1500,v3=1500,v4=1500,amp=[1,-1],f0=20,ricker=true,exponent=1,sinusoidal=false)	

	sinusoidalA = 10.*dt
	sinusoidalB = 10.*pi/(ox1 + dx1*nx1)
	sinusoidalC = 0.
	
	if length(tau1) != length(tau2) || length(tau1) != length(tau3) || length(tau1) != length(tau4)
		tau2 = 0*tau1
		tau3 = 0*tau1
		tau4 = 0*tau1
	end
	if length(v1) != length(v2) || length(v1) != length(v3) || length(v1) != length(tau4)
		v2 = 99999*v1
		v3 = 99999*v1
		v4 = 99999*v1
	end
	if length(amp) != length(v1)
		amp = 1 + 0*v1
	end
	if length(f0) != length(v1)
		f0 = 20 + 0*v1
	end

	d = zeros(nt,nx1,nx2,nx3,nx4)	
	nf = 4*nextpow2(nt) + 1
	dw = 2.*pi/nf/dt
	nw = round(Int,floor(nf/2)) + 1
	t = collect(0:1:nt-1)'*dt
	x1 = Array(Float32,1,nx1,nx2,nx3,nx4)
	x2 = Array(Float32,1,nx1,nx2,nx3,nx4)
	x3 = Array(Float32,1,nx1,nx2,nx3,nx4)
	x4 = Array(Float32,1,nx1,nx2,nx3,nx4)
	for ix1=1:nx1 
		for ix2=1:nx2 
			for ix3=1:nx3 
				for ix4=1:nx4  
					x1[1,ix1,ix2,ix3,ix4] = (ix1*dx1 + ox1)^exponent
					x2[1,ix1,ix2,ix3,ix4] = (ix2*dx2 + ox2)^exponent
					x3[1,ix1,ix2,ix3,ix4] = (ix3*dx3 + ox3)^exponent
					x4[1,ix1,ix2,ix3,ix4] = (ix4*dx4 + ox4)^exponent
				end
			end
		end
	end

	D = zeros(Complex{Float64},nf,nx1,nx2,nx3,nx4)
	for ievent=1:length(tau1)
		t,wav = SeisWavelets.Ricker(fc=f0[ievent],dt=dt)
		nwav = length(wav)
		wav = cat(2,wav',zeros(1,nf-length(wav)))
		Wav = fft(wav',1)
		delay = dt*(round(Int,nwav/2)+1)
		for iw=1:nw
			w = (iw-1)*dw
			shift = tau1[ievent] + tau2[ievent] + tau3[ievent] + tau4[ievent] + x1/v1[ievent] + x2/v2[ievent] + x3/v3[ievent] + x4/v4[ievent] - delay
			if sinusoidal
				shift += sinusoidalA*sin(sinusoidalB*x1 + sinusoidalC) + sinusoidalA*sin(sinusoidalB*x2 + sinusoidalC) + sinusoidalA*sin(sinusoidalB*x3 + sinusoidalC) + sinusoidalA*sin(sinusoidalB*x4 + sinusoidalC)
			end	
			if (ricker)
				D[iw,:,:,:,:] += amp[ievent]*Wav[iw]*exp(-1im*w*shift)
			else
				D[iw,:,:,:,:] += amp[ievent]*exp(-1im*w*(shift + delay))
			end	
		end
	end
	# symmetries
	for iw=nw+1:nf
		D[iw,:,:,:,:] = conj(D[nf-iw+2,:,:,:,:])
	end 
	d = ifft(D,1)
	d = real(d[1:nt,:,:,:,:])
	extent = Extent(convert(Int32,nt),convert(Int32,nx1),convert(Int32,nx2),convert(Int32,nx3),convert(Int32,nx4),
		   convert(Float32,ot),convert(Float32,ox1),convert(Float32,ox2),convert(Float32,ox3),convert(Float32,ox4),
		   convert(Float32,dt),convert(Float32,dx1),convert(Float32,dx2),convert(Float32,dx3),convert(Float32,dx4),
		   "Time","ix1","ix2","ix3","ix4",
		   "s","index","index","index","index",
		   "")	
	if extent.n5 == 1 && extent.n4 == 1 && extent.n3 == 1 && extent.n2 == 1 
		d = reshape(d,round(Int,extent.n1))
	elseif extent.n5 == 1 && extent.n4 == 1 && extent.n3 == 1
		d = reshape(d,round(Int,extent.n1),round(Int,extent.n2))
	elseif extent.n5 == 1 && extent.n4 == 1
		d = reshape(d,round(Int,extent.n1),round(Int,extent.n2),round(Int,extent.n3))
	elseif extent.n5 == 1
		d = reshape(d,round(Int,extent.n1),round(Int,extent.n2),round(Int,extent.n3),round(Int,extent.n4))
	else
		d = reshape(d,round(Int,extent.n1),round(Int,extent.n2),round(Int,extent.n3),round(Int,extent.n4),round(Int,extent.n5))
	end
	return d,extent
end
