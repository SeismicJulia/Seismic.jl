function ApplyModelWeights(m,param::Dict{Any,Any})

	wm = get(param,"wm",1.)
	return m.*wm

end

function ApplyModelWeights(m,wm,h1::Array{Header,1},h2::Array{Header,1},param::Dict{Any,Any})

	m_out = ApplyModelWeights(m,{"wm"=>wm})
	return m_out,h1

end

function ApplyModelWeights(m::ASCIIString,d::ASCIIString,param=Dict())

	adj = get(param,"adj",true)
	wm = get(param,"wm","wm")
	ntrace = get(param,"ntrace",100000)
	param["f"] = [ApplyModelWeights]
	param["group"] = "some"
	param["ntrace"] = ntrace
	if (adj==true)
		SeisProcess(d,wm,m,param)
	else
		SeisProcess(m,wm,d,param)
	end

end

function ApplyModelWeights(m::Array{ASCIIString,1},d::Array{ASCIIString,1},param=Dict())

	for j = 1 : length(m)
		ApplyModelWeights(m[j],d[j],param)
	end

end


