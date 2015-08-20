function CGStep(m1,m2,a,b)
# m1 = a*m1 + b*m2
	return a*m1 + b*m2
end

function CGStep(m1::ASCIIString,m2::ASCIIString,a,b)
# m1 = a*m1 + b*m2

	itrace = 1
	status = false
	tmp_filename = join(["tmp_CGStep_",string(int(rand()*100000))])
	while status == false
		ma,h,status = SeisRead(m1,"some",["tracenum"],itrace,10000)
		mb,h,status = SeisRead(m2,"some",["tracenum"],itrace,10000)
		m = a*ma + b*mb
		SeisWrite(tmp_filename,m,h,itrace)
		itrace += 10000
	end
	SeisCopy(tmp_filename,m1)
	SeisRemove(tmp_filename)
	
end

function CGMult(m::ASCIIString,m1::ASCIIString,m2::ASCIIString,a,b)
# m = a*m1*m2 + b

	itrace = 1
	status = false
	tmp_filename = join(["tmp_CGMult_",string(int(rand()*100000))])
	while status == false
		ma,h,status = SeisRead(m1,"some",["tracenum"],itrace,10000)
		mb,h,status = SeisRead(m2,"some",["tracenum"],itrace,10000)
		mc = a*ma.*mb + b
		SeisWrite(tmp_filename,mc,h,itrace)
		itrace += 10000
	end
	SeisCopy(tmp_filename,m)
	SeisRemove(tmp_filename)
	
end

function CGSparseNorm(P::ASCIIString,m::ASCIIString)
# P = abs(m)/maximum(m[:])

	maxval = SeisMaximum(m)
	itrace = 1
	status = false
	while status == false
		ma,h,status = SeisRead(m,"some",["tracenum"],itrace,10000)
		ma = abs(ma)/maxval
		SeisWrite(P,abs(ma)/maxval,h,itrace)
		itrace += 10000
	end
	
end

function SeisMaximum(m::ASCIIString)

	maxval = 0.
	status = false
	itrace = 1
	while status == false
		ma,hm,status = SeisRead(m,"some",["tracenum"],itrace,10000)
		if (maximum(abs(ma[:])) > maxval) 
			maxval = maximum(abs(ma[:]))
		end
		itrace += 10000
	end
	return maxval
end

