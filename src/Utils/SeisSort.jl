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

	filename_h = success(`grep "headers=" $filename`) ? chomp(readall(`grep "headers=" $filename` |> `tail -1` |> `awk '{print substr($1,10,length($1)-10) }' `)) : "NULL"
	stream_h = open(filename_h)
	seek(stream_h, header_count["n1"])
	nt = read(stream_h,Int32)
	nx = int32(filesize(stream_h)/(4*length(names(Header))))
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
	FetchHeaders(in,out,p,int32(nx),ntrace)
	Seismic.FetchTraces(in,out,ntrace)
	param["f"]=[UpdateHeader]
	tmp = join(["tmp_SeisSort_",string(int(rand()*100000))])
	SeisProcessHeaders(out,tmp,param)
	cp(tmp_h,out_h); rm(tmp_h);

end

function FetchHeaders(in::ASCIIString,out::ASCIIString,p::Array{Int32,1},nx::Int32;ntrace=1000)

	h = Header[]
	itrace = 1
	for j = 1:nx
		append!(h,SeisReadHeaders(in,group="some",itrace=p[j],ntrace=1))
		if (length(h) == ntrace)
			SeisWriteHeaders(out,h,itrace=itrace)
			itrace += ntrace
			h = Header[]
		end
	end 
	if (length(h) > 0) 
		SeisWriteHeaders(out,h,itrace=itrace)
	end	
	
end
