function SeisToSegy(in,out;su=true)
    
    if (su==true)
	file_hsize = 0
    else
	file_hsize = 900
	# add commands here to read text and binary headers out.thead and
        # out.bhead and write them to out
    end
    
    filename_headers = ParseHeaderName(in)
    extent = ReadTextHeader(in)
    n1 = extent.n1
    nx = extent.n2*extent.n3*extent.n4*extent.n5
    
    stream = open(out,"w")
    total = 60 + n1
    h_segy = Array(SegyHeader,1)
    h_seis = Array(Header,1)
    h1 = Header[]
    push!(h1,Seismic.InitSeisHeader())
    h1[1].o1 = extent.o1
    h1[1].d1 = extent.d1
    h1[1].n1 = extent.n1
    for j = 1 : nx
	if filename_headers != "NULL"
	    d,h1,e = SeisRead(in,itrace=j,ntrace=1,group="some")
	else
	    println("j=",j)
	    d,e = SeisRead(in,itrace=j,ntrace=1,group="some")
	    println("size(d)=",size(d))	
	    h1[1].tracenum = j
	end	
	h_segy[1] = MapHeaders(h1,j,"SeisToSegy")     
	position = file_hsize + total*(j-1)*4 + segy_count["trace"]
	seek(stream,position)
	write(stream,convert(Array{Float32,1},d[:]))
	PutSegyHeader(stream,h_segy[1],n1,file_hsize,j)
    end
    close(stream)
    
end
