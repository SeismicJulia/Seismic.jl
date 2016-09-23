function SeisMute(in;offset=[0.],tmute=0.,vmute=1500.,taper=0.1,dt=0.001)

	nt,nx = size(in[:,:])
	out = copy(in[:,:])
	for it = 1:nt
		for ix = 1:nx
			t = sqrt(tmute^2 + (offset[ix]/vmute).^2)
			if ((it-1)*dt <= t ) 
					if ((it-1)*dt > t - taper)
						out[it,ix] *= 1 - (t - (it-1)*dt)/taper
					else
						out[it,ix] = 0.			
					end
			end
		end
	end
	
	return out
end

function SeisMute(in,h::Array{Header,1};offset=[0],tmute=0.,vmute=1500.,taper=0.1,dt=0.001)

	nt,nx = size(in)
	offset = Float32[]
	for itrace = 1 : nx
		#push!(offset,sqrt((h[itrace].gx - h[itrace].sx)^2 + (h[itrace].gy - h[itrace].sy)^2 ))
		push!(offset,abs(h[itrace].h))
	end
	out = SeisMute(in,offset=offset,tmute=tmute,vmute=vmute,taper=taper,dt=h[1].d1)	
	
	return out,h
end

function SeisMute(m::String,d::String,adj;tmute=0.,vmute=1500.,taper=0.1)

	@compat parameters = Dict(:tmute=>tmute,:vmute=>vmute,:taper=>taper)
	if (adj==true)
		SeisProcess(d,m,[SeisMute],[parameters];key=["imx"])
	else
		SeisProcess(m,d,[SeisMute],[parameters];key=["imx"])
	end

end
