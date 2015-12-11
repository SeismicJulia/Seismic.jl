function ApplyDataWeights(d,param::Dict{Any,Any})

	wd = get(param,"wd",1.)
	return d.*wd

end

function ApplyDataWeights(d,wd,h1::Array{Header,1},h2::Array{Header,1},param::Dict{Any,Any})

	d_out = ApplyDataWeights(d,{"wd"=>wd})
	return d_out,h1

end

function ApplyDataWeights(m::ASCIIString,d::ASCIIString,param=Dict())

	adj = get(param,"adj",true)
	wd = get(param,"wd","wd")
	ntrace = get(param,"ntrace",100000)
	param["f"] = [ApplyDataWeights]
	param["group"] = "some"
	param["ntrace"] = ntrace
	if (adj==true)
		SeisProcess(d,wd,m,param)
	else
		SeisProcess(m,wd,d,param)
	end

end

function ApplyDataWeights(m::Array{ASCIIString,1},d::Array{ASCIIString,1},param=Dict())

	for j = 1 : length(m)
		ApplyDataWeights(m[j],d[j],param)
	end     

end                                                                                                                     

