include("Header.jl")

function SeisRead(filename,group="all",key=["imx","imy"],itrace=1,ntrace=100)

    filename_h = join([filename ".seish"])
    stream_h = open(filename_h)
    seek(stream_h, header_count["n1"])
    nt = read(stream_h,Int32)
    nhead = 27
    curr = zeros(length(key),1)
    prev = 1*curr
    nx = int(filesize(stream_h)/124) - itrace + 1
    if (group == "all")
        ntrace = nx
    elseif (group == "gather")
        for j=itrace:itrace+nx-1
            h1 = GrabHeader(stream_h,j)
        	for ikey=1:length(key)
                    curr[ikey] = getfield(h1,symbol(key[ikey]))
        	end
        	if curr != prev && j > itrace
        	    nx = j - itrace
        	    break
        	end
        	prev = 1*curr
        end
    else
    	nx > ntrace ? nx = 1*ntrace : 0
    end
    h = Array(Header,nx)
    d = zeros(nt,nx)

    filename_d = join([filename ".seisd"])
    stream_d = open(filename_d)
    position_d = 4*nt*(itrace-1)
    seek(stream_d,position_d)
    position_h = 4*(itrace-1)
    seek(stream_h,position_h)
    for j = itrace : itrace + nx - 1
    	d[:,j-itrace+1] = read(stream_d,Float32,nt)
        h[j-itrace+1] = GrabHeader(stream_h,j)
    end
    status = eof(stream_h)
    if (nx == 0)
    	status = true
    end
    close(stream_h)
    close(stream_d)
    return d,h,status
    
end
