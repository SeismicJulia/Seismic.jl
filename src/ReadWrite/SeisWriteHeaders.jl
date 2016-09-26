function SeisWriteHeaders(filename,h;itrace=1,update_tracenum=true)

    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
    filename_d = join([DATAPATH filename "@data@"])	
    filename_h = join([DATAPATH filename "@headers@"])	
    if (itrace==1)
	stream_hout = open(filename_h,"w")
    else
	stream_hout = open(filename_h,"a")
    end
    nx = length(h)
    h1 = Header32Bits[]
    for j = itrace : itrace + nx - 1
	if update_tracenum == true
	    h[j - itrace + 1].tracenum = j
	end	
	h2 = HeaderToBits(h[j - itrace + 1])
	append!(h1,h2)
    end
    write(stream_hout,h1)
    close(stream_hout)
    
end
