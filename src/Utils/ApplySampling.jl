function ApplySampling(in,wd)

	for itrace = 1 : size(in,2)
		in[:,itrace] = in[:,itrace].*wd[itrace]
	end
	return in
	
end

function ApplySampling(in::ASCIIString,wd)

	itrace = 1
	status = false
	tmp_filename = join(["tmp_ApplySampling_",string(int(rand()*100000))])
	while status == false
		d,h,status = SeisRead(in,"some",["tracenum"],itrace,1)
		d *= wd[itrace] 
		SeisWrite(tmp_filename,d,h,itrace)
		itrace += 1
	end
	SeisCopy(tmp_filename,in)
	SeisRemove(tmp_filename)
	
end

