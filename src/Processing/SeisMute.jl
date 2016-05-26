function SeisMute(in;offset=[0.],tmute=0.,vmute=1500.,taper=0.1,dt=0.001)

	nt,nx = size(in)
	out = copy(in)
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
		push!(offset,sqrt((h[itrace].gx - h[itrace].sx)^2 + (h[itrace].gy - h[itrace].sy)^2 ))
	end
	out = SeisMute(in,offset=offset,tmute=tmute,vmute=vmute,taper=taper,dt=h[1].dt)	
	
	return out,h
end

function SeisMute(m::ASCIIString,d::ASCIIString,param=Dict())

	adj = get(param,"adj",true)
	wd = get(param,"wd","wd")
	ntrace = get(param,"ntrace",100000)
	param["f"] = [SeisMute]
	param["group"] = "some"
	param["ntrace"] = ntrace
	if (adj==true)
		SeisProcess(d,m,param)
	else
		SeisProcess(m,d,param)
	end

end

function SeisMute(m::Array{ASCIIString,1},d::Array{ASCIIString,1},param=Dict())

	for j = 1 : length(m)
		SeisMute(m[j],d[j],param)
	end     

end                                                                                                                     
