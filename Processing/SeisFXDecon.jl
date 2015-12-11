function SeisFXDecon(in,param=Dict())

	dt = get(param,"dt",0.002)
	filter_length = get(param,"filter_length",10)
	mu = get(param,"mu",0.01)
	nt = size(in,1)
	nx = size(in,2)
	nf = nt
	dw = 2.*pi/nf/dt
	nw = int(nf/2) + 1
	fmax = get(param,"fmax",int(floor(0.5/dt)))
	if(fmax*dt*nf < nw) 
		iw_max = int(floor(fmax*dt*nf))
	else 
		iw_max = int(floor(0.5/dt))
	end
	D = fft(in,1)
	Df = zeros(Complex64,size(D))
	Db = zeros(Complex64,size(D))
	for iw=1:iw_max
		y = vec(D[iw,:])
		yf = fxdecon(y,"f",filter_length,mu);
		yb = fxdecon(y,"b",filter_length,mu);
		Df[iw,:] = yf
		Db[iw,:] = yb
	end
	# symmetries
	for iw=nw+1:nf
		Df[iw,:] = conj(Df[nf-iw+2,:])
		Db[iw,:] = conj(Db[nf-iw+2,:])
	end 
	d = ( ifft(Df,1) + ifft(Db,1) )/ 2
	out = real(d[1:nt,:])

	return out
end

function fxdecon(d,direction,filter_length,mu)
	nx = length(d);
	if (direction=="b")
		y  = d[1:nx-filter_length];
		C  = d[2:nx-filter_length+1];
		R  = d[nx-filter_length+1:nx];
		M = hankel(C,R);
	else
		y  = d[filter_length+1:nx,1];
		C  = d[filter_length:nx-1];
		R = flipud(d[1:filter_length]);
		M = toeplitz(C,R);
	end	
	B = M'*M;  
	beta = B[1,1]*mu/100
	ab = (B + beta*eye(filter_length))\M'*y
	temp = M*ab
	if (direction=="b")
		d = [temp ; zeros(filter_length)]
	else
		d = [zeros(filter_length) ; temp]
	end	

	return d
end

function hankel(c, r)

	nc = length(r)
	nr = length(c)
	if (nc > 1)
		c = [ c ; r[2:nc] ]
	end
	m = c[ones(nr) * [1:nc]' + [0:nr-1]*ones(nc)']
	m = reshape(m,nr,nc)

	return m
end

function toeplitz(c, r)

	nc = length(r)
	nr = length(c)
	m = zeros(Complex64,nr, nc);
	for i = 1:min(nc, nr)
		m[i:nr,i] = c[1:nr-i+1]
	end
	for i = 1:min(nr, nc-1)
		m[i,i+1:nc] = r[2:nc-i+1];
	end

	return m
end
