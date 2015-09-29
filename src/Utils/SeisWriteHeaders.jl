include("Header.jl")

function SeisWriteHeaders(filename,h,param=Dict())

	itrace = get(param,"itrace",1)

	filename_h = join([filename ".seish"])
	if (itrace==1)
		stream_hout = open(filename_h,"w")
	else
		stream_hout = open(filename_h,"a")
	end
	nx = length(h)
	h1 = Header32Bits[]
	for j = itrace : itrace + nx - 1
		h2 = HeaderToBits(h[j - itrace + 1])
		append!(h1,h2)
	end
	write(stream_hout,h1)
	close(stream_hout)
end
