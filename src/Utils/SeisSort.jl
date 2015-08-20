include("Header.jl")

function SeisSort(in, out, key=["imx","imy"], rev=false)

    filename_h = join([in ".seish"])
    stream_h = open(filename_h)
    seek(stream_h, header_count["n1"])
    nt = read(stream_h,Int32)
    nhead = 27
    nx = int(filesize(stream_h)/124)
    h = Header[]    
    # find min and max for each key
    h1 = GrabHeader(stream_h,1)
    minval = Array(Float32,length(key))
    
    for ikey=1:length(key)
    	minval[ikey] = getfield(h1,symbol(key[ikey]))
    end
    for j=2:nx
    	h1 = GrabHeader(stream_h,j)
        for ikey=1:length(key)
        	key_val = abs(getfield(h1,symbol(key[ikey])))
            if (key_val < minval[ikey])
        		minval[ikey] = key_val
            end
        end
    end
    mykey = vec(zeros(Float64,nx))
    seekstart(stream_h)
    for j=1:nx
    	h1 = GrabHeader(stream_h,j)
        for ikey=1:length(key)
        	mykey[j] += getfield(h1,symbol(key[ikey]))*(10^(6*(length(key)-ikey)))
        end
    end
	close(stream_h)
	p = sortperm(mykey,rev=rev)    
	for j = 1:nx
		d,h1,status = SeisRead(in,"some",["tracenum"],p[j],1)
		h1[1].tracenum = convert(Int32,j)
		SeisWrite(out,d,h1,j)
    end   
    
end
