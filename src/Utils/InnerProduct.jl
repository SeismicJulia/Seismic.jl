function InnerProduct(in1,in2)
	
	return convert(Float32,real(sum(conj(in1[:]).*in2[:])))
	
end

function InnerProduct(in1::ASCIIString,in2::ASCIIString;ntrace=10000)
	# Inner Product of two vectors with same length
	
		ip = 0.0
		itrace_in = 1
		nx = GetNumTraces(in1)
		while itrace_in <= nx
			d1,h1 = SeisRead(in1,group="some",itrace=itrace_in,ntrace=ntrace)
			d2,h2 = SeisRead(in2,group="some",itrace=itrace_in,ntrace=ntrace)
			ip += InnerProduct(d1,d2)		
			num_traces_in = size(d1,2)
			itrace_in += num_traces_in
		end

		return ip
	
end

function InnerProduct(in1::Array{ASCIIString,1},in2::Array{ASCIIString,1};ntrace=10000)

		ip = float32(0)
		for j = 1 : length(in1)
			ip += InnerProduct(in1[j],in2[j],ntrace=ntrace)
		end

		return ip
	
end

