"""
**SeisWrite**

*Write seismic data in .seis format*

**IN**

* filename
* d: data as 2d array
* h: headers as 1d array with elements of type Header
* itrace=1

**OUT**

*Credits: Aaron Stanton, 2015*

"""

include("Header.jl")

function SeisWrite(filename,d,h;itrace=1)

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
	h1 = Header32Bits[]
	for j = itrace : itrace + nx - 1
		h[j - itrace + 1].tracenum = j 
		h2 = HeaderToBits(h[j - itrace + 1])
		append!(h1,h2)
	end
	write(stream_hout,h1)
	close(stream_hout)
end
