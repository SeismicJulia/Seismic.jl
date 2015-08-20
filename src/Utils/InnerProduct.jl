function InnerProduct(in1,in2)
# Inner Product of two vectors with same length
	return convert(Float64,real(sum(conj(in1[:]).*in2[:])))
end

function InnerProduct(in1::ASCIIString,in2::ASCIIString)
# Inner Product of two vectors with same length
	itrace = 1
	status = false
	ip = 0.
	while status == false
		d1,h,status = SeisRead(in1,"some",["tracenum"],itrace,1)
		d2,h,status = SeisRead(in1,"some",["tracenum"],itrace,1)
		ip += sum(conj(d1[:]).*d2[:])
		itrace += 1
	end
	return convert(Float64,real(ip))
end
