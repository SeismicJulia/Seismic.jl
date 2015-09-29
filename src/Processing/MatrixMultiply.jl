function MatrixMultiply(in,param::Dict{Any,Any})

	adj = get(param,"adj",false)
	matrix = get(param,"matrix",eye(length(in)))
	if (adj)
		out = matrix'*in
	else
		out = matrix*in
	end

	return out
end

function MatrixMultiply(in,h::Array{Header,1},param=Dict())

	adj = get(param,"adj",false)
	matrix = get(param,"matrix",eye(length(in)))
	
	if (adj)
		out = zeros(Float32,size(matrix,2))
	else
		out = zeros(Float32,size(matrix,1))
	end
	for j = 1 : length(h)
		if (adj)
			out[:,j] = matrix'*in[:,j]
			h[j].n1 = size(matrix,2)
		else
			out[:,j] = matrix*in[:,j]
			h[j].n1 = size(matrix,1)
		end
	end
	return out,h
end

function MatrixMultiply(m::ASCIIString,d::ASCIIString,param=Dict())

	adj = get(param,"adj",false)
	ntrace = get(param,"ntrace",100000)
	param["f"] = [MatrixMultiply]
	param["group"] = "some"
	param["ntrace"] = ntrace
	if (adj==true)
		SeisProcess(d,m,param)
	else
		SeisProcess(m,d,param)
	end

end
