function CalculateSampling(in,h=[],param=Dict())

	cutoff = get(param,"cutoff",1e-10)
	itrace = 1
	wd = zeros(Float32,size(in))
	for itrace = 1 : size(in[:,:],2)
		a = sum(in[:,itrace].*in[:,itrace])
		if (a > cutoff) 
			wd[:,itrace] = 1.
		end
	end
	return wd,h;
end

function CalculateSampling(in::ASCIIString,wd::ASCIIString,param=Dict())
	# calculate sampling operator (1's for live traces, 0's for missing traces)
	
	ntrace = get(param,"ntrace",100000)
	cutoff = get(param,"cutoff",1e-10)
	param["f"] = [CalculateSampling]
	param["group"] = "some"
	param["ntrace"] = ntrace
	param["cutoff"] = cutoff
	SeisProcess(in,wd,param)		
		
end

