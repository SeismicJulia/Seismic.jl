function SeisProcess(in,out,f,param,group="some",key=["imx","imy"],ntrace=100)	
# Run processing flows that read and write from disk
#
# group="all"  will process all the traces from the input at once
# 		"some" will process ntrace traces from the input at once
# 		"gather" will use the "key" parameter and the headers to process gathers
#        if using group="gather" then "key" can be set to a vector of headers to 
#        key on. 
#
# f is a function that has the following syntax: d2,h2 = f(d1,h1,param), where 
#
# param is a dictionary (Dict) of parameters for the function.
# note that f can be a vector of functions. They will be executed sequentially on
# the same group of traces.
#

	if (group=="all")
		d1,h1,status = SeisRead(in,group,key,1,ntrace)
    	for ifunc = 1 : length(f)
    		func = f[ifunc]
    		d2,h2 = func(d1,h1,param)
    		d1 = copy(d2)
    		h1 = copy(h2)
    	end
		SeisWrite(out,d1,h1)
    else
        itrace_in = 1
        itrace_out = 1
        status = false
    	while status == false
    		d1,h1,status = SeisRead(in,group,key,itrace_in,ntrace)
    		num_traces_in = size(d1,2)
    		for ifunc = 1 : length(f)
    			func = f[ifunc]
    		    d2,h2 = func(d1,h1,param)
    			d1 = copy(d2)
    			h1 = copy(h2)
    		end
    		num_traces_out = size(d1,2)
    		for itrace = 1 : num_traces_out
    			h1[itrace].tracenum = itrace_out + itrace - 1
    		end
		    SeisWrite(out,d1,h1,itrace_out)
		    itrace_in += num_traces_in
		    itrace_out += num_traces_out
        end
    end
    
end
