function CGStep(m1,m2,h1=Header[],h2=Header[],param=Dict())
	# m1 = a*m1 + b*m2
	a = get(param,"a",1.)
	b = get(param,"b",1.)
	return a*m1 + b*m2,h1
end

function CGStep(m1::ASCIIString,m2::ASCIIString,param=Dict())
	# m1 = a*m1 + b*m2	

	a = get(param,"a",1.)
	b = get(param,"b",1.)
	param["group"]="some"
	param["f"]=[CGStep]
	ntrace = get(param,"ntrace",100000)
	tmp = join(["tmp_CGStep_",string(int(rand()*100000))])
	SeisProcess(m1,m2,tmp,param)
	SeisCopy(tmp,m1)
	SeisRemove(tmp)
end

function CGStep(m1::Array{ASCIIString,1},m2::Array{ASCIIString,1},param=Dict())
	# m1 = a*m1 + b*m2	
	a = get(param,"a",[1. 1. 1.])
	b = get(param,"b",[1. 1. 1.])
	for j = 1 : length(m1)
		CGStep(m1[j],m2[j],{"a"=>a[j],"b"=>b[j]})
	end
	
end

function CGSparseNorm(m,h,param=Dict())

	maxval = get(param,"maxval",1.)
	P = abs(m)/maxval
	return P,h

end

function CGSparseNorm(m::ASCIIString,P::ASCIIString,param=Dict())
	# P = abs(m)/maximum(m[:])

	param["maxval"] = SeisMaximum(m,param)
	ntrace = get(param,"ntrace",100000)
	param["f"] = [CGSparseNorm]
	param["group"] = "some"
	param["ntrace"] = ntrace
	SeisProcess(m,P,param)

end

function SeisMaximum(m::ASCIIString,param=Dict())

	ntrace = get(param,"ntrace",100000)
	param["ntrace"] = ntrace

	maxval = 0.
	itrace_in = 1
	nx = GetNumTraces(m)
	while itrace_in <= nx
		m1,h1 = SeisRead(m,["group"=>"some","itrace"=>itrace_in,"ntrace"=>ntrace])
		if (maximum(abs(m1[:])) > maxval) 
			maxval = maximum(abs(m1[:]))
		end
		num_traces_in = size(m1,2)
		itrace_in += num_traces_in
	end
	return maxval
end

