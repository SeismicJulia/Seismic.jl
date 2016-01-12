include("Header.jl")

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
* extent: extent of the data (try _names(Extent)_ to see the information this contains)

*Credits: AS, 2015*

"""
function SeisRead(filename;group="all",key=["imx","imy"],itrace=1,ntrace=10000)

	filename_data = chomp(readall(`grep "in=" $filename` |> `tail -1` |> `awk '{print substr($1,5,length($1)-5) }' `))
	filename_headers = success(`grep "headers=" $filename`) ? chomp(readall(`grep "headers=" $filename` |> `tail -1` |> `awk '{print substr($1,10,length($1)-10) }' `)) : "NULL"
	extent = ReadTextHeader(filename)
	stream_d = open(filename_data)
	dtype = chomp(readall(`grep "data_format" $filename` |> `tail -1` |> `awk '{print substr($1,14,length($1)-14)}'` ))
	dtype = dtype == "native_float" ? Float32 : Complex{Float32}
	esize = int(chomp(readall(`grep "esize" $filename` |> `tail -1` |> `awk '{print substr($1,7,length($1))}'` )))
	total = int(filesize(stream_d)/esize)
	nhead = length(names(Header))
	curr = zeros(length(key),1)
	prev = 1*curr
	nx = extent.n2*extent.n3*extent.n4*extent.n5
	if filename_headers != "NULL"
		stream_h = open(filename_headers)
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
			ntrace = nx > ntrace ? ntrace : nx
		else
			ntrace = nx > ntrace ? ntrace : nx
		end
		position_h = 4*nhead*(itrace-1)
		seek(stream_h,position_h)
		h1 = read(stream_h,Header32Bits,nhead*ntrace)
		h1 = reshape(h1,nhead,int(ntrace))
		h = Header[]
		for j = 1 : ntrace
			h = push!(h,BitsToHeader(h1[:,j]))    	
		end    
		close(stream_h)
	else
		ntrace = nx
	end
	stream_d = open(filename_data)
	position_d = 4*extent.n1*(itrace-1)
	seek(stream_d,position_d)
	d = read(stream_d,dtype,extent.n1*ntrace)
	if extent.n5 == 1 && extent.n4 == 1 && extent.n3 == 1 && extent.n2 == 1 
		d = reshape(d,int(extent.n1))
	elseif extent.n5 == 1 && extent.n4 == 1 && extent.n3 == 1
		d = reshape(d,int(extent.n1),int(extent.n2))
	elseif extent.n5 == 1 && extent.n4 == 1
		d = reshape(d,int(extent.n1),int(extent.n2),int(extent.n3))
	elseif extent.n5 == 1
		d = reshape(d,int(extent.n1),int(extent.n2),int(extent.n3),int(extent.n4))
	else
		d = reshape(d,int(extent.n1),int(extent.n2),int(extent.n3),int(extent.n4),int(extent.n5))
	end
	close(stream_d)
	if filename_headers != "NULL"
		return d,h,extent
	else
		return d,extent	
	end
end
