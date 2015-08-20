include("SegyStruct.jl")
include("Header.jl")

function SeisToSegy(in,out,su=true)

	if (su==true)
		file_hsize = 0
	else
		file_hsize = 900
		# add commands here to read text and binary headers out.thead and out.bhead and write them to out
	end
    stream_h = open(join([in ".seish"]))
    seek(stream_h, file_hsize + header_count["n1"])
    nt = read(stream_h,Int32)
    nx = int(filesize(stream_h)/124)
    close(stream_h)
    
    stream = open(out,"w")
    total = 60 + nt
    d = zeros(Float32,nt,1)
    h_segy = Array(SegyHeader,1)
    h_seis = Array(Header,1)
    for j = 1 : nx
    	d,h1,status = SeisRead(in,"some",["tracenum"],j,1)
        h_segy[1] = MapHeaders(h1,j,"SeisToSegy")     
    	position = file_hsize + total*(j-1)*4 + segy_count["trace"]
    	seek(stream,position)
    	write(stream,convert(Array{Float32,1},d[:]))
        PutSegyHeader(stream,h_segy[1],nt,file_hsize,j)
    end
    close(stream)
	
end


