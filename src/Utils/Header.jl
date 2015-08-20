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

header_count = Dict{String,Int32}()
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
    position = 124*(j-1)
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
    position = 124*(j-1)
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

