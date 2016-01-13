function SeisProcess(in::ASCIIString,out::ASCIIString,operators,parameters;group="some",key=[],ntrace=10000)	
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
		d1,h1,e1 = SeisRead(in,group=group,key=key,itrace=1,ntrace=ntrace)
		for j = 1 : length(operators)
			op = operators[j]
			d2 = op(d1,h1,parameters[j])
			d1 = copy(d2)
		end
		SeisWrite(out,d1,h1,e1)
	else
		itrace_in = 1
		itrace_out = 1
		nx = GetNumTraces(in)
		while itrace_in <= nx
			d1,h1,e1 = SeisRead(in,group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			num_traces_in = size(d1,2)
			for j = 1 : length(operators)
				op = operators[j]
				d2 = op(d1,h1,parameters[j])
				d1 = copy(d2)
			end
			num_traces_out = size(d1,2)
			SeisWrite(out,d1,h1,e1,itrace=itrace_out)
			itrace_in += num_traces_in
			itrace_out += num_traces_out
		end
	end

end

function SeisProcess(in::Array{ASCIIString,1},out::Array{ASCIIString,1},operators,parameters;group="some",key=[],ntrace=10000)
 	
	for j = 1 : length(in)
		SeisProcess(in[j],out[j],param)
	end

end

function SeisProcess(in1::ASCIIString,in2::ASCIIString,out::ASCIIString,operators,parameters;group="some",key=[],ntrace=10000)
	f = get(param,"f",[fft_op])
	group = get(param,"group","some")
	key = get(param,"key",["imx","imy"])
	ntrace = get(param,"ntrace",100)

	if (group=="all")
		d1,h1 = SeisRead(in1,group=group,key=key,itrace=1,ntrace=ntrace)
		d2,h2 = SeisRead(in2,group=group,key=key,itrace=1,ntrace=ntrace)
		for ifunc = 1 : length(f)
			func = f[ifunc]
			d3,h3 = func(d1,d2,h1,h2,param)
			d1 = copy(d3)
			h1 = copy(h3)
		end
		SeisWrite(out,d1,h1)
	else
		itrace_in = 1
		itrace_out = 1
		nx = GetNumTraces(in1)
		while itrace_in <= nx
			d1,h1 = SeisRead(in1,group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			d2,h2 = SeisRead(in2,group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			num_traces_in = size(d1,2)
			for ifunc = 1 : length(f)
				func = f[ifunc]
				d3,h3 = func(d1,d2,h1,h2,param)
				d1 = copy(d3)
				h1 = copy(h3)
			end
			num_traces_out = size(d1,2)
			SeisWrite(out,d1,h1,itrace=itrace_out)
			itrace_in += num_traces_in
			itrace_out += num_traces_out
		end
	end

end
