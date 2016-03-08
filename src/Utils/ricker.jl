function ricker(f0,dt)
	nw=2.2/f0/dt
	nw=2*int(floor(nw/2))+1
	nc=int(nw/2)
	w = zeros(nw)
	k=[1:1:nw]
  k=vec(k)
	alpha = (nc-k+1)*f0*dt*pi
	beta=alpha.^2
	w = (1.-beta.*2).*exp(-beta)
	return w
end
