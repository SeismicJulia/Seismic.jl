include("Header.jl")

function SeisWrite(filename,d,h,itrace=1)
    
    filename_d = join([filename ".seisd"])
    filename_h = join([filename ".seish"])
    if (itrace==1)
    	stream_dout = open(filename_d,"w")
    	stream_hout = open(filename_h,"w")
    else
    	stream_dout = open(filename_d,"a")
    	stream_hout = open(filename_h,"a")
    end
    write(stream_dout,convert(Array{Float32,1},d[:]))
    close(stream_dout)
    nx = size(d,2)
    for j = itrace : itrace + nx - 1
    	h[j - itrace + 1].tracenum = j
    	PutHeader(stream_hout,h[j - itrace + 1],j)
    end
    close(stream_hout)
end
