function SeisNMO(in,h,param=Dict())

	nt,nx = size(in)
	dt = h[1].d1

	offset = Float32[]
	for itrace = 1 : nx
		push!(offset,h[itrace].h)
	end

	tnmo = get(param,"tnmo",[0:5:nt-1]*dt)
	vnmo = get(param,"vnmo",tnmo*0 + 1500)
	max_stretch = get(param,"max_stretch",100)
	# interpolate tau,v pairs to match sampling of input data 
	if (length(vnmo) == 1)
		tnmo = convert(Float64,tnmo);
		vnmo = convert(Float64,vnmo);
		ti = [0:1:nt-1]*dt
		vi = ones(1,nt)*vnmo
	else
		tnmo = convert(Array{Float64,1},vec(tnmo))
		vnmo = convert(Array{Float64,1},vec(vnmo))
		ti = [0:1:nt-1]*dt
		g = InterpIrregular(tnmo, vnmo, BCnan, InterpLinear)
		vi = g[ti]
	end
	out = zeros(size(in))
	M = zeros(nt,1)
	for it = 1:nt
		for ix = 1:nx
			time = sqrt(ti[it]^2 + (offset[ix]/vi[it]).^2)
			stretch = (time-ti[it])/(ti[it]+1e-10)
			if (stretch<max_stretch/100)
				its = time/dt+1
				it1 = floor(time/dt+1)
				it2 = it1+1
				a = its-it1
				if (it2 <= nt) 
					out[it,ix] = (1-a)*in[it1,ix]+a*in[it2,ix]
				end
			end
		end
	end
	return out,h
end
