"""
**SegyToSeis**

*Convert SEGY or SU data to seis format.*

**IN**

* `filename_in::AbstractString`
* `filename_out::AbstractString`
* format="segy" (or "su")
* swap_bytes=true (flag to swap bytes)
* input_type="ibm" (or "ieee")

**OUT**

*Credits: AS, 2015*

"""

function SegyToSeis(filename_in::AbstractString,filename_out::AbstractString;format="segy",swap_bytes=true,input_type="ibm")

	if (format=="su")
		file_hsize = 0
	else
		file_hsize = 3600
		# add commands here to read text and binary headers and write them out to
		# filename_out.thead and filename_out.bhead
		stream     = open(filename_in)
		position = 3200
		seek(stream, position)
		fh = GrabFileHeader(stream)
		ntfh = swap_bytes == true ? bswap(fh.netfh) : fh.netfh
		if ntfh == -1
			error("add instructions to deal with variable extended text header")
		end
		if ntfh == 0
			file_hsize = 3600
		elseif ntfh > 0
			# file_hsize = 3200 * (ntfh+1) + 400
			file_hsize = 3200 * 1 + 400
		else
			error("unknown data format")
		end
	end
	stream = open(filename_in)
	seek(stream, segy_count["ns"] + file_hsize)
	if (swap_bytes==true)
		nt = bswap(read(stream,Int16))
	else
		nt = read(stream,Int16)
	end
	total = 60 + nt
	nx = round(Int,(filesize(stream)-file_hsize)/4/total)
	println("number of traces: ",nx)
	println("number of samples per trace: ",nt)
	h_segy = Array{SegyHeader}(1)
	h_seis = Array{Header}(1)
	seek(stream,file_hsize + segy_count["trace"])
	h_segy[1] = GrabSegyHeader(stream,swap_bytes,nt,file_hsize,1)
	dt = h_segy[1].dt/1000000
	extent = Extent(convert(Int32,nt),convert(Int32,nx),convert(Int32,1),convert(Int32,1),convert(Int32,1),
		   convert(Float32,0),convert(Float32,1),convert(Float32,0),convert(Float32,0),convert(Float32,0),
		   convert(Float32,dt),convert(Float32,1),convert(Float32,1),convert(Float32,1),convert(Float32,1),
		   "Time","Trace Number","","","",
		   "s","index","","","",
		   "")
	for j=1:nx
		position = file_hsize + total*(j-1)*4 + segy_count["trace"]
		seek(stream,position)
		if (input_type == "ieee")
			d = read(stream,Float32,nt)
		else
			d = read(stream,IBMFloat32,nt)
		end
		if (swap_bytes==true && input_type == "ieee")
			d = bswap_vector(d)
		end
		if (input_type != "ieee")
			d = convert(Array{Float32,1},d)
		end
		h_segy[1] = GrabSegyHeader(stream,swap_bytes,nt,file_hsize,j)
		h_seis[1] = MapHeaders(h_segy,j,"SegyToSeis")
		SeisWrite(filename_out,d,h_seis,extent,itrace=j)
	end
	close(stream)

end

function bswap_vector(a)
	for i = 1 : length(a)
		a[i] = bswap(a[i]);
	end
	return a
end
