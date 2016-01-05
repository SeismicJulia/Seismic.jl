"""
**SeisRead**

*Read seismic data in .seis format*

**IN**   

* filename
* group="all" ("some" or "gather")
* key=["imx","imy"]
* itrace=1
* ntrace=10000

**OUT**  

* d: data as 2d array 
* h: headers as 1d array 

*Credits: Aaron Stanton, 2015*

"""

include("Header.jl")

function SeisRead(filename;group="all",key=["imx","imy"],itrace=1,ntrace=10000)

	filename_h = join([filename ".seish"])
	stream_h = open(filename_h)
	seek(stream_h, header_count["n1"])
	nt = read(stream_h,Int32)
	nhead = length(names(Header))
	curr = zeros(length(key),1)
	prev = 1*curr
	nx = int(filesize(stream_h)/(4*length(names(Header)))) - itrace + 1
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
		ntrace = nx > ntrace ? ntrace : nx
	else
		ntrace = nx > ntrace ? ntrace : nx
	end
	
	filename_d = join([filename ".seisd"])
	stream_d = open(filename_d)
	position_d = 4*nt*(itrace-1)
	seek(stream_d,position_d)
	position_h = 4*nhead*(itrace-1)
	seek(stream_h,position_h)

	d = read(stream_d,Float32,nt*ntrace)
	h1 = read(stream_h,Header32Bits,nhead*ntrace)
	d = reshape(d,int(nt),int(ntrace))
	h1 = reshape(h1,nhead,int(ntrace))
	h = Header[]
	for itrace = 1 : ntrace
		h = push!(h,BitsToHeader(h1[:,itrace]))    	
	end    

	close(stream_h)
	close(stream_d)
	return d,h

end
