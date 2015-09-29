include("Header.jl")

function SeisWindow(in,out,param=Dict())

	key = get(param,"key",[])
	minval = get(param,"minval",[])
	maxval = get(param,"maxval",minval)
	SeisWindowHeaders(in,out,param)
	filename = join([in ".seish"])
	stream = open(join([in ".seish"]))
	seek(stream, header_count["o1"])
	ot = read(stream,Float32)
	seek(stream, header_count["n1"])
	nt = read(stream,Int32)
	seek(stream, header_count["d1"])
	dt = read(stream,Float32)
	close(stream)
	tmin = ot
	tmax = ot + dt*(nt-1)
	key2 = Any[]
	minval2 = Any[]
	maxval2 = Any[]
	for ikey=1:length(key)
		if key[ikey] == "t"
			tmin = minval[ikey]
			tmax = maxval[ikey]
		end
	end
	itmin = int32((tmin - ot)/dt + 1)
	if (itmin < 1)
		itmin = 1
	end
	itmax = int32((tmax - ot)/dt + 1)  
	if (itmax) > nt
		itmax = nt
	end
	param["itmin"] = itmin
	param["itmax"] = itmax
	FetchTraces(in,out,param);
	param["f"]=[UpdateHeader]
	tmp = join(["tmp_SeisWindow_",string(int(rand()*100000))])
	SeisProcessHeaders(out,tmp,param)
	cp(join([tmp ".seish"]),join([out ".seish"])); rm(join([tmp ".seish"]));
end

function FetchTraces(in,out,param=Dict())

	ntrace = get(param,"ntrace",500)
	NX = GetNumTraces(out)
	itrace = int32(1)
	stream_in  = open(join([in ".seisd"]))
	stream_out = open(join([out ".seisd"]),"w")

	while itrace <= NX
		if (itrace > 1)
			stream_out = open(join([out ".seisd"]),"a")
		end
		nx = NX - itrace + 1
		ntrace = nx > ntrace ? ntrace : nx
		h = SeisReadHeaders(out,["group"=>"some","itrace"=>itrace,"ntrace"=>ntrace])
		nt = h[1].n1
		itmin = get(param,"itmin",int32(1))
		itmax = get(param,"itmax",int32(nt))
		d = zeros(Float32,itmax-itmin+1,ntrace)
		SeekTraces!(d,stream_in,h,itmin,itmax,nt,int32(ntrace))
		write(stream_out,d)
		close(stream_out)
		itrace += ntrace
	end
	close(stream_in)
end	

function SeekTraces!{T}(d::AbstractArray{T,2},stream_in::IOStream,h::Array{Header,1},itmin::Int32,itmax::Int32,nt::Int32,ntrace::Int32)

	d1 = zeros(Float32,nt)
	for ix = 1 : ntrace
			position = 4*nt*(h[ix].tracenum-1)
			seek(stream_in,position)
			d1 = read(stream_in,Float32,nt)
			for it = itmin : itmax
				d[it - itmin + 1,ix] = d1[it]	
			end
	end
	nothing
end

function UpdateHeader(h,param=Dict())

	itrace = get(param,"itrace",1)
	itmin = get(param,"itmin",1)
	itmax = get(param,"itmax",h[1].n1)
	ot = (itmin-1)*h[1].d1 + h[1].o1
	nt = itmax - itmin + 1
	for j = 1 : length(h)
		h[j].o1 = ot
		h[j].n1 = nt
		h[j].tracenum = j + itrace - 1
	end

	return h
end
