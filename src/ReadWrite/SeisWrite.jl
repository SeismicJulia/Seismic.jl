"""
**SeisWrite**

*Write seismic data in .seis format*

**IN**

* filename
* d: data
* h: headers as 1d array with elements of type Header
* extent: extent of the data (try _names(Extent)_ to see the information this contains)
* itrace=1

**OUT**

*Credits: AS, 2015*

"""
function SeisWrite(filename,d,h::Array{Header,1},extent::Extent;itrace=1)

    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
    filename_d = join([DATAPATH filename "@data@"])
    filename_h = join([DATAPATH filename "@headers@"])	
    if (itrace==1)
	WriteTextHeader(filename,extent,"native_float",4,filename_d,filename_h)
	stream_dout = open(filename_d,"w")
	stream_hout = open(filename_h,"w")
    else
	stream_dout = open(filename_d,"a")
	stream_hout = open(filename_h,"a")
    end
    write(stream_dout,convert(Array{Float32,1},d[:]))
    close(stream_dout)
    nt=size(d,1)
    nx = size(reshape(d,nt,:),2)
    h1 = Header32Bits[]
    for j = itrace : itrace + nx - 1
	h[j - itrace + 1].tracenum = j 
	h2 = HeaderToBits(h[j - itrace + 1])
	append!(h1,h2)
    end
    write(stream_hout,h1)
    close(stream_hout)
end

function SeisWrite(filename,d,extent::Extent)
    
    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
    filename_d = join([DATAPATH filename "@data@"])
    WriteTextHeader(filename,extent,"native_float",4,filename_d,"NULL")
    stream_dout = open(filename_d,"w")
    write(stream_dout,convert(Array{Float32,1},d[:]))
    close(stream_dout)
end
