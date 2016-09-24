#include("../ReadWrite/Header.jl")

"""
**SeisWindow**
*Window a seis file using header words.*
**IN**   
* in
* out
* key=[]
* minval=[]
* maxval=[]
* note that windowing along the time axis is achieved by using the key "t".
**OUT**  
*Credits: AS, 2015*
"""

function SeisWindow(in::String,out::String;key=[],minval=[],maxval=[])
    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
	extent = ReadTextHeader(in)
	tmin = extent.o1
	tmax = extent.o1 + extent.d1*(extent.n1-1)
	for ikey=1:length(key)
		if key[ikey] == "t"
			tmin = minval[ikey]
			tmax = maxval[ikey]
		end
	end
	itmin = convert(Int32,round((tmin - extent.o1)/extent.d1) + 1)
	if (itmin < 1)
		itmin = 1
	end
	itmax = convert(Int32,round((tmax - extent.o1)/extent.d1) + 1)  
	if (itmax) > extent.n1
		itmax = extent.n1
	end
	SeisWindowHeaders(in,out;key=key,minval=minval,maxval=maxval,tmin=tmin,tmax=tmax)
	FetchTraces(in,out;itmin=itmin,itmax=itmax)
	tmp = join(["tmp_SeisWindow_",string(round(Int,rand()*100000))])
	@compat SeisProcessHeaders(out,tmp,[UpdateHeader],[Dict(:itmin=>itmin,:itmax=>itmax)])
	filename_h_tmp = join([DATAPATH tmp "@headers@"])	
	filename_h_out = join([DATAPATH out "@headers@"])	
	cp(filename_h_tmp,filename_h_out,remove_destination=true);
	rm(filename_h_tmp);
	#rm(tmp);
end

function FetchTraces(in::String,out::String;ntrace=500,itmin=round(Int,1),itmax=round(Int,9999999999))

	NX = GetNumTraces(out)
	itrace = round(Int,1)
	filename_data_in = ParseDataName(in)
    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
	filename_data_out = join([DATAPATH out "@data@"])
	stream_in  = open(filename_data_in)
	stream_out  = open(filename_data_out,"w")
	extent = ReadTextHeader(in)
	nt = extent.n1
	if itmin < 1
		itmin = round(Int,1)
	end
	if itmax > nt
		itmax = round(Int,nt)
	end
	while itrace <= NX
		if (itrace > 1)
			stream_out = open(filename_data_out,"a")
		end
		nx = NX - itrace + 1
		ntrace = nx > ntrace ? ntrace : nx
		h = SeisReadHeaders(out,group="some",itrace=itrace,ntrace=ntrace)
		d = zeros(Float32,itmax-itmin+1,ntrace)
		SeekTraces!(d,stream_in,h,itmin,itmax,nt,round(Int,ntrace))
		write(stream_out,d)
		close(stream_out)
		itrace += ntrace
	end
	close(stream_in)
end	

function SeekTraces!{T}(d::AbstractArray{T,2},stream_in::IOStream,h::Array{Header,1},itmin,itmax,nt,ntrace)

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

function UpdateHeader(h;itrace=1,itmin=1,itmax=9e9)
	ot = (itmin-1)*h[1].d1 + h[1].o1
	nt = itmax - itmin + 1
	for j = 1 : length(h)
		h[j].o1 = ot
		h[j].n1 = nt
		h[j].tracenum = j + itrace - 1
	end

	return h
end
