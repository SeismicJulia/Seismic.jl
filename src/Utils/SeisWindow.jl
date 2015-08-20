include("Header.jl")

function SeisWindow(in,out,key=["imx"],minval=[-9999999],maxval=minval)

	filename_h = join([in ".seish"])
    stream_h = open(filename_h)
    seek(stream_h, header_count["n1"])
    nt = read(stream_h,Int32)
    seek(stream_h, header_count["d1"])
    dt = read(stream_h,Float32)
    
    tmin = 0
    tmax = 999999999
    key2 = Any[]
    minval2 = Any[]
    maxval2 = Any[]
    for ikey=1:length(key)
    	if key[ikey] != "t"
    		key2 = push!(key2,key[ikey])
    		minval2 = push!(minval2,minval[ikey])
    		maxval2 = push!(maxval2,maxval[ikey])
    	else
    		tmin = minval[ikey]
    		tmax = maxval[ikey]
    	end
    end
    itmin = int(floor(tmin/dt) + 1)
    if (itmin < 1)
    	itmin = 1
    end
    itmax = int(floor(tmax/dt) + 1)  
    if (itmax) > nt
    	itmax = nt
    end
    
    nhead = 27
    nx = int(filesize(stream_h)/124)
    h = Header[]
    for j=1:nx
    	h1 = GrabHeader(stream_h,j)
    	keep = true
        for ikey=1:length(key2)
        	key_val = getfield(h1,symbol(key2[ikey]))
            if (key_val < minval2[ikey] || key_val > maxval2[ikey])
        		keep = false
            end
        end
		if (keep==true)
			h = push!(h,h1)
		end
    end
    close(stream_h)
    nx = length(h)
	if (nx==0)
		error(join(["no traces selected from ",in]))
	end

    filename_d = join([in ".seisd"])
    stream_d = open(filename_d)
    d = zeros(nt,1)
	#println("nx=",nx)
	h1 = Array(Header,1)
    for j=1:nx
    	position_d = 4*nt*(h[j].tracenum - 1)
    	seek(stream_d,position_d)
    	d[:] = read(stream_d,Float32,nt)
    	d1 = d[itmin:itmax]
    	h1[1] = h[j]
    	h1[1].tracenum = convert(Int32,j)
    	h1[1].o1 = convert(Float32,tmin)
    	h1[1].n1 = convert(Int32,itmax-itmin+1)
		SeisWrite(out,d1,h1,j)
    end
    close(stream_d)
    
end
