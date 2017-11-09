function SeisSemblance(in,param=Dict())

	dt = get(param,"dt",0.002)
	dx = get(param,"dx",10)
	nt,nx = size(in)
	offset = get(param,"offset",collect(0:1:nx-1)*dx)
	vmin = get(param,"vmin",1500)
	vmax = get(param,"vmax",3500)
	nv = get(param,"nv",nx)
	v = linspace(vmin,vmax,nv)
	tau = collect(0:4:nt-1)*dt
	ntau = length(tau)
	L = get(param,"L",5)

	S = zeros(ntau,nv)
	for it = 1:ntau
		for iv = 1:nv
			t = sqrt( tau[it].^2 + (offset[iv]/v[iv]).^2 )
			s = zeros(2*L+1,nx)
			for ig = -L:L
				ts = t + (ig-1)*dt
				for ix = 1:nx
					is = ts/dt+1
					i1 = convert(Int64,floor(is))
					i2 = i1 + 1
					if (i1>=1 && i2<=nt)
						a = is-i1
						s[ig+L+1,ix] = (1.-a)*in[i1,ix] + a*in[i2,ix]   #Grab sample with Linear interpolation
					end
				end
			end
			#s = s.*H
			s1  = sum( (sum(s,2)).^2)
			s2  = sum( sum(s.^2))
			S[it,iv] = abs(s1-s2)
		end
	end
	S = S/maximum(S[:])

	return S,tau,v
end
