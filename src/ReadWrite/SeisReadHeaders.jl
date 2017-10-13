"""
**SeisReadHeaders**
*Read seismic headers in seis format (binary file). 

**
**IN**   
* filename
* group="all" ("some" or "gather")
* key=["imx","imy"]
* itrace=1
* ntrace=10000
**OUT**  
* h: headers as 1d array 
*Credits: AS, 2015*
"""

function SeisReadHeaders(filename;group="all",key=[],itrace=1,ntrace=100)

	filename_h = ParseHeaderName(filename)
	stream_h = open(filename_h)
	@compat nhead = length(fieldnames(Header))
	curr = zeros(length(key),1)
	prev = 1*curr
	@compat nx = round(Int,filesize(stream_h)/(4*length(fieldnames(Header)))) - itrace + 1
	if (group == "all")
		ntrace = nx
	elseif (group == "gather")
		for j=itrace:itrace+nx-1
			h1 = GrabHeader(stream_h,j)
			for ikey=1:length(key)
				curr[ikey] = getfield(h1,Symbol(key[ikey]))
			end
			if curr != prev && j > itrace
				nx = j - itrace
				break
			end
			prev = 1*curr
		end
	else
		ntrace = nx > ntrace ? ntrace : nx
	end

	position_h = 4*nhead*(itrace-1)
	seek(stream_h,position_h)

	h1 = read(stream_h,Header32Bits,nhead*ntrace)
	h1 = reshape(h1,nhead,round(Int,ntrace))
	h = Header[]
	for itrace = 1 : ntrace
		h = push!(h,BitsToHeader(h1[:,itrace]))    	
	end    

	close(stream_h)
	return h

end
