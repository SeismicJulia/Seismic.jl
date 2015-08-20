function CalculateSampling(in)

    itrace = 1
    wd = zeros(Float32,size(in))
    for itrace = 1 : size(in,2)
    	a = sum(in[:,itrace].*in[:,itrace])
    	if (a > 0.00001) 
    		wd[:,itrace] = 1.
    	end
    end
	return wd;
end

function CalculateSampling(in::ASCIIString,wd::ASCIIString)
# calculate sampling operator (1's for live traces, 0's for missing traces)
    itrace = 1
    status = false
    while status == false
    	d,h,status = SeisRead(in,"some",["tracenum"],itrace,10000)
	wd1 = d;
	#println("itrace="size(d))
	#println(size(wd))
	for ix = 1 : size(d,2)
		a = sum(d[:,ix].*d[:,ix])
    		if (a < 0.00001) 
    			wd1[itrace,:] = 0.
    		else
    			wd1[itrace,:] = 1.
		end	
	end	
	SeisWrite(wd,wd1,h,itrace)
	itrace += 10000
    end
end

