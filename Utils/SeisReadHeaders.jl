include("Header.jl")

function SeisReadHeaders(filename,param=Dict())

	group = get(param,"group","all")
	key = get(param,"key",["imx","imy"])
	itrace = get(param,"itrace",1)
	ntrace = get(param,"ntrace",100)

	filename_h = join([filename ".seish"])
	stream_h = open(filename_h)
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
	else
		ntrace = nx > ntrace ? ntrace : nx
	end

	position_h = 4*nhead*(itrace-1)
	seek(stream_h,position_h)

	h1 = read(stream_h,Header32Bits,nhead*ntrace)
	h1 = reshape(h1,nhead,int(ntrace))
	h = Header[]
	for itrace = 1 : ntrace
		h = push!(h,BitsToHeader(h1[:,itrace]))    	
	end    

	close(stream_h)
	return h

end
