type Header
	tracenum::Int32
	o1::Float32
	n1::Int32
	d1::Float32
	sx::Float32
	sy::Float32
	gx::Float32
	gy::Float32
	mx::Float32
	my::Float32
	hx::Float32
	hy::Float32
	h::Float32
	az::Float32
	ang::Float32
	isx::Int32
	isy::Int32
	igx::Int32
	igy::Int32
	imx::Int32
	imy::Int32
	ihx::Int32
	ihy::Int32
	ih::Int32
	iaz::Int32
	iang::Int32
	selev::Float32
	gelev::Float32
	sstat::Float32
	gstat::Float32
	trid::Int32
end

header_count = Dict{AbstractString,Int32}()
header_count["tracenum"] = 0
header_count["o1"]    = 4
header_count["n1"]    = 8
header_count["d1"]    = 12
header_count["sx"]    = 16
header_count["sy"]    = 20
header_count["gx"]    = 24
header_count["gy"]    = 28
header_count["mx"]    = 32
header_count["my"]    = 36
header_count["hx"]    = 40
header_count["hy"]    = 44
header_count["h"]     = 48
header_count["az"]    = 52
header_count["ang"]   = 56
header_count["isx"]   = 60
header_count["isy"]   = 64
header_count["igx"]   = 68
header_count["igy"]   = 72
header_count["imx"]   = 76
header_count["imy"]   = 80
header_count["ihx"]   = 84
header_count["ihy"]   = 88
header_count["ih"]    = 92
header_count["iaz"]   = 96
header_count["iang"]  = 100
header_count["selev"] = 104
header_count["gelev"] = 108
header_count["sstat"] = 112
header_count["gstat"] = 116
header_count["trid"]  = 120

function InitSeisHeader()
	h = Header(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
	0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
	0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
	0.0)
	return h
end

function GrabHeader(stream,j)
	position = 4*length(names(Header))*(j-1)
	seek(stream,position)
	h = InitSeisHeader()
	h.tracenum = read(stream,typeof(h.tracenum))
	h.o1    = read(stream,typeof(h.o1))
	h.n1    = read(stream,typeof(h.n1))
	h.d1    = read(stream,typeof(h.d1)) 
	h.sx    = read(stream,typeof(h.sx))
	h.sy    = read(stream,typeof(h.sy))
	h.gx    = read(stream,typeof(h.gx))
	h.gy    = read(stream,typeof(h.gy))
	h.mx    = read(stream,typeof(h.mx))
	h.my    = read(stream,typeof(h.my))
	h.hx    = read(stream,typeof(h.hx))
	h.hy    = read(stream,typeof(h.hy))
	h.h     = read(stream,typeof(h.h))
	h.az    = read(stream,typeof(h.az))
	h.ang   = read(stream,typeof(h.ang))
	h.isx   = read(stream,typeof(h.isx))
	h.isy   = read(stream,typeof(h.isy))
	h.igx   = read(stream,typeof(h.igx))
	h.igy   = read(stream,typeof(h.igy))
	h.imx   = read(stream,typeof(h.imx))
	h.imy   = read(stream,typeof(h.imy))
	h.ihx   = read(stream,typeof(h.ihx))
	h.ihy   = read(stream,typeof(h.ihy))
	h.ih    = read(stream,typeof(h.ih))
	h.iaz   = read(stream,typeof(h.iaz))
	h.iang   = read(stream,typeof(h.iang))
	h.selev = read(stream,typeof(h.selev))
	h.gelev = read(stream,typeof(h.gelev))
	h.sstat = read(stream,typeof(h.sstat))
	h.gstat = read(stream,typeof(h.gstat))
	h.trid  = read(stream,typeof(h.trid))
	return h
end

function PutHeader(stream,h,j)
	position = 4*length(names(Header))*(j-1)
	seek(stream,position)
	write(stream,h.tracenum)
	write(stream,h.o1)
	write(stream,h.n1)
	write(stream,h.d1) 
	write(stream,h.sx)
	write(stream,h.sy)
	write(stream,h.gx)
	write(stream,h.gy)
	write(stream,h.mx)
	write(stream,h.my)
	write(stream,h.hx)
	write(stream,h.hy)
	write(stream,h.h)
	write(stream,h.az)
	write(stream,h.ang)
	write(stream,h.isx)
	write(stream,h.isy)
	write(stream,h.igx)
	write(stream,h.igy)
	write(stream,h.imx)
	write(stream,h.imy)
	write(stream,h.ihx)
	write(stream,h.ihy)
	write(stream,h.ih)
	write(stream,h.iaz)
	write(stream,h.iang)
	write(stream,h.selev)
	write(stream,h.gelev)
	write(stream,h.sstat)
	write(stream,h.gstat)
	write(stream,h.trid)
end

bitstype 32 Header32Bits

function BitsToHeader(h_in)
	h = InitSeisHeader()
	h.tracenum = reinterpret(typeof(h.tracenum),h_in[1])
	h.o1    = reinterpret(typeof(h.o1),h_in[2])
	h.n1    = reinterpret(typeof(h.n1),h_in[3])
	h.d1    = reinterpret(typeof(h.d1),h_in[4]) 
	h.sx    = reinterpret(typeof(h.sx),h_in[5])
	h.sy    = reinterpret(typeof(h.sy),h_in[6])
	h.gx    = reinterpret(typeof(h.gx),h_in[7])
	h.gy    = reinterpret(typeof(h.gy),h_in[8])
	h.mx    = reinterpret(typeof(h.mx),h_in[9])
	h.my    = reinterpret(typeof(h.my),h_in[10])
	h.hx    = reinterpret(typeof(h.hx),h_in[11])
	h.hy    = reinterpret(typeof(h.hy),h_in[12])
	h.h     = reinterpret(typeof(h.h),h_in[13])
	h.az    = reinterpret(typeof(h.az),h_in[14])
	h.ang   = reinterpret(typeof(h.ang),h_in[15])
	h.isx   = reinterpret(typeof(h.isx),h_in[16])
	h.isy   = reinterpret(typeof(h.isy),h_in[17])
	h.igx   = reinterpret(typeof(h.igx),h_in[18])
	h.igy   = reinterpret(typeof(h.igy),h_in[19])
	h.imx   = reinterpret(typeof(h.imx),h_in[20])
	h.imy   = reinterpret(typeof(h.imy),h_in[21])
	h.ihx   = reinterpret(typeof(h.ihx),h_in[22])
	h.ihy   = reinterpret(typeof(h.ihy),h_in[23])
	h.ih    = reinterpret(typeof(h.ih),h_in[24])
	h.iaz   = reinterpret(typeof(h.iaz),h_in[25])
	h.iang   = reinterpret(typeof(h.iang),h_in[26])
	h.selev = reinterpret(typeof(h.selev),h_in[27])
	h.gelev = reinterpret(typeof(h.gelev),h_in[28])
	h.sstat = reinterpret(typeof(h.sstat),h_in[29])
	h.gstat = reinterpret(typeof(h.gstat),h_in[30])
	h.trid  = reinterpret(typeof(h.trid),h_in[31])
	return h
end

function HeaderToBits(h);
	h_out = [reinterpret(Header32Bits,h.tracenum);
	reinterpret(Header32Bits,h.o1);
	reinterpret(Header32Bits,h.n1);
	reinterpret(Header32Bits,h.d1);
	reinterpret(Header32Bits,h.sx);
	reinterpret(Header32Bits,h.sy);
	reinterpret(Header32Bits,h.gx);
	reinterpret(Header32Bits,h.gy);
	reinterpret(Header32Bits,h.mx);
	reinterpret(Header32Bits,h.my);
	reinterpret(Header32Bits,h.hx);
	reinterpret(Header32Bits,h.hy);
	reinterpret(Header32Bits,h.h);
	reinterpret(Header32Bits,h.az);
	reinterpret(Header32Bits,h.ang);
	reinterpret(Header32Bits,h.isx);
	reinterpret(Header32Bits,h.isy);
	reinterpret(Header32Bits,h.igx);
	reinterpret(Header32Bits,h.igy);
	reinterpret(Header32Bits,h.imx);
	reinterpret(Header32Bits,h.imy);
	reinterpret(Header32Bits,h.ihx);
	reinterpret(Header32Bits,h.ihy);
	reinterpret(Header32Bits,h.ih);
	reinterpret(Header32Bits,h.iaz);
	reinterpret(Header32Bits,h.iang);
	reinterpret(Header32Bits,h.selev);
	reinterpret(Header32Bits,h.gelev);
	reinterpret(Header32Bits,h.sstat);
	reinterpret(Header32Bits,h.gstat);
	reinterpret(Header32Bits,h.trid)]
	return h_out
end

function GetNumTraces(in)

	filename_h = success(`grep "headers=" $in`) ? chomp(readall(`grep "headers=" $in` |> `tail -1` |> `awk '{print substr($1,10,length($1)-10) }' `)) : "NULL"
	nhead = length(names(Header))
	stream_h = open(filename_h)
	nx = int(filesize(stream_h)/(nhead*4))
	close(stream_h)
	return nx

end        

function ExtractHeader(h::Array{Header,1},key::ASCIIString)
	
	keytype = eval(parse("typeof(Seismic.InitSeisHeader().$(key))"))
	out = keytype[]
	for ix = 1 : length(h)
	push!(out,getfield(h[ix],symbol(key)))
	end
	return out
	
end

type Extent
        n1::Int32
        n2::Int32
        n3::Int32
        n4::Int32
        n5::Int32
        o1::Float32
        o2::Float32
        o3::Float32
        o4::Float32
        o5::Float32
        d1::Float32
        d2::Float32
        d3::Float32
        d4::Float32
        d5::Float32
        label1::ASCIIString
        label2::ASCIIString
        label3::ASCIIString
        label4::ASCIIString
        label5::ASCIIString
        unit1::ASCIIString
        unit2::ASCIIString
        unit3::ASCIIString
        unit4::ASCIIString
        unit5::ASCIIString
        title::ASCIIString
end

function ReadTextHeader(filename)
	n1 = success(`grep "n1" $filename`) ? int(chomp(readall(`grep "n1" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	n2 = success(`grep "n2" $filename`) ? int(chomp(readall(`grep "n2" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	n3 = success(`grep "n3" $filename`) ? int(chomp(readall(`grep "n3" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	n4 = success(`grep "n4" $filename`) ? int(chomp(readall(`grep "n4" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	n5 = success(`grep "n5" $filename`) ? int(chomp(readall(`grep "n5" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	o1 = success(`grep "o1" $filename`) ? float(chomp(readall(`grep "o1" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	o2 = success(`grep "o2" $filename`) ? float(chomp(readall(`grep "o2" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	o3 = success(`grep "o3" $filename`) ? float(chomp(readall(`grep "o3" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	o4 = success(`grep "o4" $filename`) ? float(chomp(readall(`grep "o4" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	o5 = success(`grep "o5" $filename`) ? float(chomp(readall(`grep "o5" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	d1 = success(`grep "d1" $filename`) ? float(chomp(readall(`grep "d1" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	d2 = success(`grep "d2" $filename`) ? float(chomp(readall(`grep "d2" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	d3 = success(`grep "d3" $filename`) ? float(chomp(readall(`grep "d3" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	d4 = success(`grep "d4" $filename`) ? float(chomp(readall(`grep "d4" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	d5 = success(`grep "d5" $filename`) ? float(chomp(readall(`grep "d5" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	label1 = success(`grep "label1" $filename`) ? chomp(readall(`grep "label1" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	label2 = success(`grep "label2" $filename`) ? chomp(readall(`grep "label2" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	label3 = success(`grep "label3" $filename`) ? chomp(readall(`grep "label3" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	label4 = success(`grep "label4" $filename`) ? chomp(readall(`grep "label4" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	label5 = success(`grep "label5" $filename`) ? chomp(readall(`grep "label5" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	unit1 = success(`grep "unit1" $filename`) ? chomp(readall(`grep "unit1" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	unit2 = success(`grep "unit2" $filename`) ? chomp(readall(`grep "unit2" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	unit3 = success(`grep "unit3" $filename`) ? chomp(readall(`grep "unit3" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	unit4 = success(`grep "unit4" $filename`) ? chomp(readall(`grep "unit4" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	unit5 = success(`grep "unit5" $filename`) ? chomp(readall(`grep "unit5" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	title = success(`grep "title" $filename`) ? chomp(readall(`grep "title" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	extent = Extent(convert(Int32,n1),convert(Int32,n2),convert(Int32,n3),convert(Int32,n4),convert(Int32,n5),
		   convert(Float32,o1),convert(Float32,o2),convert(Float32,o3),convert(Float32,o4),convert(Float32,o5),
		   convert(Float32,d1),convert(Float32,d2),convert(Float32,d3),convert(Float32,d4),convert(Float32,d5),
		   label1,label2,label3,label4,label5,
		   unit1,unit2,unit3,unit4,unit5,
		   title)	
	return extent
end

function WriteTextHeader(filename,extent,format,esize,filename_d,filename_h)
	# write the text header
	stream = open(filename, "w")	
	write(stream,join(["	n1=",extent.n1,"\n"]))
	write(stream,join(["	n2=",extent.n2,"\n"]))
	write(stream,join(["	n3=",extent.n3,"\n"]))
	write(stream,join(["	n4=",extent.n4,"\n"]))
	write(stream,join(["	n5=",extent.n5,"\n"]))
	write(stream,join(["	o1=",extent.o1,"\n"]))
	write(stream,join(["	o2=",extent.o2,"\n"]))
	write(stream,join(["	o3=",extent.o3,"\n"]))
	write(stream,join(["	o4=",extent.o4,"\n"]))
	write(stream,join(["	o5=",extent.o5,"\n"]))
	write(stream,join(["	d1=",extent.d1,"\n"]))
	write(stream,join(["	d2=",extent.d2,"\n"]))
	write(stream,join(["	d3=",extent.d3,"\n"]))
	write(stream,join(["	d4=",extent.d4,"\n"]))
	write(stream,join(["	d5=",extent.d5,"\n"]))
	write(stream,join(["	label1=\"",extent.label1,"\"\n"]))
	write(stream,join(["	label2=\"",extent.label2,"\"\n"]))
	write(stream,join(["	label3=\"",extent.label3,"\"\n"]))
	write(stream,join(["	label4=\"",extent.label4,"\"\n"]))
	write(stream,join(["	label5=\"",extent.label5,"\"\n"]))
	write(stream,join(["	unit1=\"",extent.unit1,"\"\n"]))
	write(stream,join(["	unit2=\"",extent.unit2,"\"\n"]))
	write(stream,join(["	unit3=\"",extent.unit3,"\"\n"]))
	write(stream,join(["	unit4=\"",extent.unit4,"\"\n"]))
	write(stream,join(["	unit5=\"",extent.unit5,"\"\n"]))
	write(stream,join(["	title=\"",extent.title,"\"\n"]))
	write(stream,join(["	data_format=\"",format,"\"\n"]))
	write(stream,join(["	esize=",esize,"\n"]))
	write(stream,join(["	in=\"",filename_d,"\"\n"]))
	write(stream,join(["	headers=\"",filename_h,"\"\n"]))
	close(stream)
end
































