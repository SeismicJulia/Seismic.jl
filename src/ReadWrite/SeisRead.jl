"""
**SeisRead**
*Read seismic data in seis format. The format is comprised of three elements:*
* a text file with information about data extent, data and header file names
* a binary file containing data
* a binary file containing headers
**
**IN**   
* filename
* group="all" ("some" or "gather")
* key=["imx","imy"]
* itrace=1
* ntrace=10000
**OUT**  
* d: data as 2d array 
* h: headers as 1d array 
* extent: extent of the data (try _fieldnames(Extent)_ to see the information this contains)
*Credits: AS, 2015*
"""
function SeisRead(filename;group="all",key=["imx","imy"],itrace=1,ntrace=10000)

    filename_data = ParseDataName(filename)
    filename_headers = ParseHeaderName(filename)
    extent = ReadTextHeader(filename)
    #println(filename_data)
    stream_d = open(filename_data)
    dtype = ParseDataFormat(filename)
    dtype = dtype == "native_float" ? Float32 : Complex{Float32}
    esize = ParseDataESize(filename)
    total = convert(Int,filesize(stream_d)/esize)
    close(stream_d)
    @compat nhead = length(fieldnames(Header))
    curr = zeros(length(key),1)
    prev = 1*curr
    nx = extent.n2*extent.n3*extent.n4*extent.n5
    if filename_headers != "NULL"
	stream_h = open(filename_headers)
	if (group == "all")
	    ntrace = nx
	elseif (group == "gather")
	    j = 1
	    for j=itrace:nx
		h1 = GrabHeader(stream_h,j)
		for ikey=1:length(key)
		    curr[ikey] = getfield(h1,Symbol(key[ikey]))
				end
		if curr != prev && j > itrace
		    ntrace = j - itrace
		    break
		end
		prev = 1*curr
	    end
	    ntrace = j < nx ? j - itrace : j - itrace + 1
	else
	    ntrace = nx - itrace + 1 > ntrace ? ntrace : nx - itrace + 1
	end
	position_h = 4*nhead*(itrace-1)
	seek(stream_h,position_h)
	h1 = read(stream_h,Header32Bits,nhead*ntrace)
	h1 = reshape(h1,nhead,convert(Int64,ntrace))
	h = Header[]
	for j = 1 : ntrace
	    h = push!(h,BitsToHeader(h1[:,j]))    	
	end    
	close(stream_h)
    else
	ntrace = nx > ntrace ? ntrace : nx
    end
	stream_d = open(filename_data)
    position_d = 4*extent.n1*(itrace-1)
    seek(stream_d,position_d)
    d = read(stream_d,dtype,extent.n1*ntrace)
    if group=="all"
	if extent.n5 == 1 && extent.n4 == 1 && extent.n3 == 1 && extent.n2 == 1 
	    d = reshape(d,Int(extent.n1))
	elseif extent.n5 == 1 && extent.n4 == 1 && extent.n3 == 1
	    d = reshape(d,convert(Int,extent.n1),convert(Int,extent.n2))
	elseif extent.n5 == 1 && extent.n4 == 1
	    d = reshape(d,convert(Int,extent.n1),convert(Int,extent.n2),convert(Int,extent.n3))
	elseif extent.n5 == 1
	    d = reshape(d,convert(Int,extent.n1),convert(Int,extent.n2),convert(Int,extent.n3),convert(Int,extent.n4))
	else
	    d = reshape(d,convert(Int,extent.n1),convert(Int,extent.n2),convert(Int,extent.n3),convert(Int,extent.n4),convert(Int,extent.n5))
	end
    else
	d = reshape(d,convert(Int64,extent.n1),convert(Int64,ntrace))
	end
    close(stream_d)
    if filename_headers != "NULL"
	return d,h,extent
    else
	return d,extent	
    end
end
