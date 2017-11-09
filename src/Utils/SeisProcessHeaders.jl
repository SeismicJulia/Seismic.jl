function SeisProcessHeaders(in, out, functions, parameters;
                            group="some", key=[], ntrace=1000000,
                            update_tracenum=true)


    if (group=="all")
	h1 = SeisReadHeaders(in,group=group,key=key,itrace=1,ntrace=ntrace)
	for ifunc = 1 : length(functions)
	    f = functions[ifunc]
	    p = parameters[ifunc]
	    h2 = f(h1;p...)
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
	    h1 = SeisReadHeaders(in, group=group, key=key, itrace=itrace_in,
                                 ntrace=ntrace)
	    num_traces_in = length(h1)

	    for ifunc = 1 : length(functions)
		f = functions[ifunc]
		p = parameters[ifunc]

		h2 = f(h1;p...)
		h1 = copy(h2)
	    end
	    num_traces_out = length(h1)

	    SeisWriteHeaders(out, h1, itrace=itrace_out,
                             update_tracenum=update_tracenum)
	    itrace_in += num_traces_in
	    itrace_out += num_traces_out
	end
    end
end
