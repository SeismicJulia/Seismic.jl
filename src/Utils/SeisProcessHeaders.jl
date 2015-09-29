function SeisProcessHeaders(in,out,param=Dict())	
	f = get(param,"f",[])
	group = get(param,"group","some")
	key = get(param,"key",["imx","imy"])
	ntrace = get(param,"ntrace",100)

	if (group=="all")
		h1 = SeisReadHeaders(in,["group"=>group,"key"=>key,"itrace"=>1,"ntrace"=>ntrace])
		for ifunc = 1 : length(f)
			func = f[ifunc]
			h2 = func(h1,param)
			h1 = copy(h2)
		end
		SeisWriteHeaders(out,h1)
	else
		itrace_in = 1
		itrace_out = 1
		NX = GetNumTraces(in)
		while itrace_in <= NX
			nx = NX - itrace_in + 1
			ntrace = nx > ntrace ? ntrace : nx
			h1 = SeisReadHeaders(in,["group"=>group,"key"=>key,"itrace"=>itrace_in,"ntrace"=>ntrace])
			num_traces_in = length(h1)
			param["itrace"]=itrace_in
			for ifunc = 1 : length(f)
				func = f[ifunc]
				h2 = func(h1,param)
				h1 = copy(h2)
			end
			num_traces_out = length(h1)
			SeisWriteHeaders(out,h1,["itrace"=>itrace_out])
			itrace_in += num_traces_in
			itrace_out += num_traces_out
		end
	end

end
