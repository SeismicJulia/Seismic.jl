include("Header.jl")

function SeisWindowHeaders(in,out;key=[],minval=[],maxval=[],tmin=0,tmax=99999,ntrace=500)
	SeisProcessHeaders(in,out,[WindowHeaders],[[:key=>key,:minval=>minval,:maxval=>maxval]],group="some",key=key,ntrace=ntrace,update_tracenum=false)
	DATAPATH = get(ENV,"DATAPATH","./")
	filename_d_out = join([DATAPATH out "@data@"])
	filename_h_out = join([DATAPATH out "@headers@"])	
	nhead = length(names(Header))
	stream_h = open(filename_h_out)
	nx = int(filesize(stream_h)/(nhead*4))
	h = GrabHeader(stream_h,1)
	close(stream_h)
	nt = convert(Int64,round((tmax - tmin)/h.d1)) + 1
	if nt > h.n1
		nt = h.n1
	end
	extent = ReadTextHeader(in)
	extent.n1 = nt
	extent.n2 = nx
	extent.o1 = tmin
	WriteTextHeader(out,extent,"native_float",4,filename_d_out,filename_h_out)

end

function WindowHeaders(h_in;key=[],minval=[],maxval=[])

	minval = convert(Array{Float32,1},vec(minval))
	maxval = convert(Array{Float32,1},vec(maxval))
	nx = length(h_in)
	key2 = ASCIIString[]
	minval2 = Float32[]
	maxval2 = Float32[]
	for ikey=1:length(key)
		if key[ikey] != "t"
			key2 = push!(key2,key[ikey])
			minval2 = push!(minval2,minval[ikey])
			maxval2 = push!(maxval2,maxval[ikey])
		end
	end
	
	return RejectHeaders(h_in,key2,minval2,maxval2,int32(length(key2)),int32(length(h_in)))
end

function RejectHeaders(h_in::Array{Header,1},key::Array{ASCIIString,1},minval::Array{Float32,1},maxval::Array{Float32,1},nkeys::Int32,nx::Int32)
	h_out = Header[]	
	keep = true
	key_val = 0f0
	for j=1:nx
		keep = true
		for ikey=1:nkeys
			key_val = convert(Float32,getfield(h_in[j],symbol(key[ikey])))
			if (key_val < minval[ikey] || key_val > maxval[ikey])
				keep = false
			end
		end
		if (keep==true)
			h_out = push!(h_out,h_in[j])
		end
	end 	
	
	return h_out
end
