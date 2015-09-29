include("Header.jl")

function SeisWindowHeaders(in,out,param=Dict())

	ntrace = get(param,"ntrace",500)
	param["f"] = [WindowingHeaders]
	param["group"] = "some"
	param["ntrace"] = ntrace
	SeisProcessHeaders(in,out,param);

end

function WindowingHeaders(h_in,param)

	key = get(param,"key",[])
	minval = convert(Array{Float32,1},vec(get(param,"minval",[])))
	maxval = convert(Array{Float32,1},vec(get(param,"maxval",minval)))
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