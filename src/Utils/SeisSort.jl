include("Header.jl")

"""
**SeisSort**

*Sort a seis file using its header words*

* a text file with information about data extent, data and header file names
* a binary file containing data
* a binary file containing headers

**

**IN**   

* in: input filename
* out: output filename
* key=["imx","imy"]
* rev=false : sort headers in decreasing order 
* ntrace=1000 : number of traces to read at a time

**OUT**  

*Credits: AS, 2015*

"""
function SeisSort(in, out;key=["imx","imy"],rev=false,ntrace=1000)

	filename_h = success(`grep "headers=" $in`) ? chomp(readall(`grep "headers=" $in` |> `tail -1` |> `awk '{print substr($1,10,length($1)-10) }' `)) : "NULL"
	stream_h = open(filename_h)
	seek(stream_h, header_count["n1"])
	nt = read(stream_h,Int32)
	nx = convert(Int64,filesize(stream_h)/(4*length(names(Header))))
	h = Header[]    
	# find min and max for each key
	h1 = GrabHeader(stream_h,1)
	minval = Array(Float32,length(key))

	for ikey=1:length(key)
		minval[ikey] = getfield(h1,symbol(key[ikey]))
	end
	for j=2:nx
		h1 = GrabHeader(stream_h,j)
		for ikey=1:length(key)
			key_val = abs(getfield(h1,symbol(key[ikey])))
			if (key_val < minval[ikey])
				minval[ikey] = key_val
			end
		end
	end
	mykey = vec(zeros(Float64,nx))
	seekstart(stream_h)
	for j=1:nx
		h1 = GrabHeader(stream_h,j)
		for ikey=1:length(key)
			mykey[j] += (getfield(h1,symbol(key[ikey])) + minval[ikey] + 1)*(10^(6*(length(key)-ikey)))
		end
	end
	close(stream_h)
	p = convert(Array{Int32,1},sortperm(mykey,rev=rev))
	FetchHeaders(in,out,p,nx)
    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
    filename_d_out = join([DATAPATH out "@data@"])
    filename_h_out = join([DATAPATH out "@headers@"])    
    nhead = length(names(Header))
    stream_h = open(filename_h_out)
    nx = int(filesize(stream_h)/(nhead*4))
    h = GrabHeader(stream_h,1)
    close(stream_h)
    extent = ReadTextHeader(in)
    extent.n2 = nx
    extent.n3 = 1
    extent.n4 = 1
    extent.n5 = 1
    WriteTextHeader(out,extent,"native_float",4,filename_d_out,filename_h_out)
	Seismic.FetchTraces(in,out)
 	tmp = join(["tmp_SeisSort_",string(int(rand()*100000))])
    @compat SeisProcessHeaders(out,tmp,[UpdateHeader],[Dict(:itmin=>1,:itmax=>nt)])
    filename_h_tmp = join([DATAPATH tmp "@headers@"])    
    filename_h_out = join([DATAPATH out "@headers@"])    
    cp(filename_h_tmp,filename_h_out);
    rm(filename_h_tmp);
	rm(tmp);
end

function FetchHeaders(in::ASCIIString,out::ASCIIString,p::Array{Int32,1},nx;ntrace=1000)

	h = Header[]
	itrace = 1
	for j = 1:nx
		append!(h,SeisReadHeaders(in,group="some",itrace=p[j],ntrace=1))
		if (length(h) == ntrace)
			SeisWriteHeaders(out,h,itrace=itrace,update_tracenum=false)
			itrace += ntrace
			h = Header[]
		end
	end 
	if (length(h) > 0) 
		SeisWriteHeaders(out,h,itrace=itrace,update_tracenum=false)
	end	
	
end
