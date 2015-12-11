include("SegyStruct.jl")
include("Header.jl")

function SegyToSeis(filename_in,filename_out,param=Dict())
	format = get(param,"format","segy")
	swap_bytes = get(param,"swap_bytes",false)
	input_type = get(param,"input_type","ieee")

	if (format=="su" || format=="rsf")
		file_hsize = 0
	else
		file_hsize = 3600
		# add commands here to read text and binary headers and write them out to 
		# filename_out.thead and filename_out.bhead
		stream     = open(filename_in)
		position = 3200
		seek(stream, position)
		fh         = GrabFileHeader(stream)
		ntfh       = fh.netfh
		if ntfh == -1
			error("need to add command to deal with variable extend text file header")
		end
		if ntfh == 0
			file_hsize = 3600
		elseif ntfh > 0
			file_hsize = 3200 * (ntfh+1) + 400
		else
			error("unknow data format")
		end
	end
	if (format=="segy" || format=="su")
		stream = open(filename_in)
		seek(stream, segy_count["ns"] + file_hsize)
		if (swap_bytes==true)
			nt = bswap(read(stream,Int16))
		else
			nt = read(stream,Int16)
		end

		total = 60 + nt
		nx = int((filesize(stream)-file_hsize)/4/total)
		println("nx=",nx)
		println("nt=",nt)

		h_segy = Array(SegyHeader,1)
		h_seis = Array(Header,1)
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
			SeisWrite(filename_out,d,h_seis,["itrace"=>j])
		end
		close(stream)
	else # read an rsf file
		rsf_filename = chomp(readall(`grep "in" $filename_in` |> `grep "@"` |> `tail -1` |> `awk '{print substr($1,5,length($1)-5) }' `))
		nt = int(chomp(readall(`grep "n1" $filename_in` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` )))
		dt = float(chomp(readall(`grep "d1" $filename_in` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` )))
		total = nt
		stream = open(rsf_filename)
		nx = int((filesize(stream)-file_hsize)/4/total)
		println("nx=",nx)
		println("nt=",nt)

		h = Header[]
		for j=1:nx
			h1 = InitSeisHeader()
			h1.n1 = nt
			h1.d1 = dt
			h1.tracenum = j
			push!(h,h1)
		end
		d = read(stream,Float32,nt*nx)
		if (swap_bytes==true)
			d = bswap_vector(d)
		end
		SeisWriteHeaders(filename_out,h)
		close(stream)
		stream = open(join([filename_out,".seisd"]),"w")
		write(stream,d)
		close(stream)
	end


end

function bswap_vector(a)
	for i = 1 : length(a)
		a[i] = bswap(a[i]);
	end
	return a
end

