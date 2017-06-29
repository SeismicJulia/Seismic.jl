type fileHeader
	jobid  :: Int32
	linnum :: Int32
	renum  :: Int32
	ntrpe  :: Int16
	natrpe :: Int16
	dt     :: Int16
	dtfr   :: Int16
	ns     :: Int16
	nsfr   :: Int16
	fmtc   :: Int16
	expf   :: Int16
	trsc   :: Int16
	vsumc  :: Int16
	sfs    :: Int16
	sfe    :: Int16
	slen   :: Int16
	styp   :: Int16
	tnumsc :: Int16
	stalens:: Int16
	stalene:: Int16
	tyta   :: Int16
	corr   :: Int16
	rgc    :: Int16
	arm    :: Int16
	unit   :: Int16
	pol    :: Int16
	vpol   :: Int16
	fvn    :: Int16
	fltf   :: Int16
	netfh  :: Int16
end

fileHeader_count = Dict{AbstractString, Int32}()
fileHeader_count["jobid"  ] = 3200
fileHeader_count["linnum" ] = 3204
fileHeader_count["renum"  ] = 3208
fileHeader_count["ntrpe"  ] = 3212
fileHeader_count["natrpe" ] = 3214
fileHeader_count["dt"     ] = 3216
fileHeader_count["dtfr"   ] = 3218
fileHeader_count["ns"     ] = 3220
fileHeader_count["nsfr"   ] = 3222
fileHeader_count["fmtc"   ] = 3224
fileHeader_count["expf"   ] = 3226
fileHeader_count["trsc"   ] = 3228
fileHeader_count["vsumc"  ] = 3230
fileHeader_count["sfs"    ] = 3232
fileHeader_count["sfe"    ] = 3234
fileHeader_count["slen"   ] = 3236
fileHeader_count["styp"   ] = 3238
fileHeader_count["tnumsc" ] = 3240
fileHeader_count["stalens"] = 3242
fileHeader_count["stalene"] = 3244
fileHeader_count["tyta"   ] = 3246
fileHeader_count["corr"   ] = 3248
fileHeader_count["rgc"    ] = 3250
fileHeader_count["arm"    ] = 3252
fileHeader_count["unit"   ] = 3254
fileHeader_count["pol"    ] = 3256
fileHeader_count["vpol"   ] = 3258
fileHeader_count["fvn"    ] = 3500
fileHeader_count["fltf"   ] = 3502
fileHeader_count["netfh"  ] = 3504

function InitFileHeader()
	fh = fileHeader(0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,)
	return fh
end

function GrabFileHeader(stream)

	position = fileHeader_count["jobid"]
	seek(stream,position)
	fh = InitFileHeader()
	fh.jobid   = read(stream,typeof(fh.jobid))
	fh.linnum  = read(stream,typeof(fh.linnum))
	fh.renum   = read(stream,typeof(fh.renum))
	fh.ntrpe   = read(stream,typeof(fh.ntrpe))
	fh.natrpe  = read(stream,typeof(fh.natrpe))
	fh.dt      = read(stream,typeof(fh.dt))
	fh.dtfr    = read(stream,typeof(fh.dtfr))
	fh.ns      = read(stream,typeof(fh.ns))
	fh.nsfr    = read(stream,typeof(fh.nsfr))
	fh.fmtc    = read(stream,typeof(fh.fmtc))
	fh.expf    = read(stream,typeof(fh.expf))
	fh.trsc    = read(stream,typeof(fh.trsc))
	fh.vsumc   = read(stream,typeof(fh.vsumc))
	fh.sfs     = read(stream,typeof(fh.sfs))
	fh.sfe     = read(stream,typeof(fh.sfe))
	fh.slen    = read(stream,typeof(fh.slen))
	fh.styp    = read(stream,typeof(fh.styp))
	fh.tnumsc  = read(stream,typeof(fh.tnumsc))
	fh.stalens = read(stream,typeof(fh.stalens))
	fh.stalene = read(stream,typeof(fh.stalene))
	fh.tyta    = read(stream,typeof(fh.tyta))
	fh.corr    = read(stream,typeof(fh.corr))
	fh.rgc     = read(stream,typeof(fh.rgc))
	fh.arm     = read(stream,typeof(fh.arm))
	fh.unit    = read(stream,typeof(fh.unit))
	fh.pol     = read(stream,typeof(fh.pol))
	fh.vpol    = read(stream,typeof(fh.vpol))
	fh.fvn     = read(stream,typeof(fh.fvn))
	fh.fltf    = read(stream,typeof(fh.fltf))
	fh.netfh   = read(stream,typeof(fh.netfh))

	return fh
end



function PutFileHeader(stream,fh)
	position = 3600
	seek(stream,position)
	write(stream,fh.jobid)
	write(stream,fh.linnum)
	write(stream,fh.renum)
	write(stream,fh.ntrpe)
	write(stream,fh.natrpe)
	write(stream,fh.dt)
	write(stream,fh.dtfr)
	write(stream,fh.ns)
	write(stream,fh.nsfr)
	write(stream,fh.fmtc)
	write(stream,fh.expf)
	write(stream,fh.trsc)
	write(stream,fh.vsumc)
	write(stream,fh.sfs)
	write(stream,fh.sfe)
	write(stream,fh.slen)
	write(stream,fh.styp)
	write(stream,fh.tnumsc)
	write(stream,fh.stalens)
	write(stream,fh.stalene)
	write(stream,fh.tyta)
	write(stream,fh.corr)
	write(stream,fh.rgc)
	write(stream,fh.arm)
	write(stream,fh.unit)
	write(stream,fh.pol)
	write(stream,fh.vpol)
	write(stream,fh.fvn)
	write(stream,fh.fltf)
	write(stream,fh.netfh)

end




type SegyHeader
	tracl::Int32
	tracr::Int32
	fldr::Int32
	tracf::Int32
	ep::Int32
	cdp::Int32
	cdpt::Int32
	trid::Int16
	nva::Int16
	nhs::Int16
	duse::Int16
	offset::Int32
	gelev::Int32
	selev::Int32
	sdepth::Int32
	gdel::Int32
	sdel::Int32
	swdep::Int32
	gwdep::Int32
	scalel::Int16
	scalco::Int16
	sx::Int32
	sy::Int32
	gx::Int32
	gy::Int32
	counit::Int16
	wevel::Int16
	swevel::Int16
	sut::Int16
	gut::Int16
	sstat::Int16
	gstat::Int16
	tstat::Int16
	laga::Int16
	lagb::Int16
	delrt::Int16
	muts::Int16
	mute::Int16
	ns::Int16
	dt::Int16
	gain::Int16
	igc::Int16
	igi::Int16
	corr::Int16
	sfs::Int16
	sfe::Int16
	slen::Int16
	styp::Int16
	stas::Int16
	stae::Int16
	tatyp::Int16
	afilf::Int16
	afils::Int16
	nofilf::Int16
	nofils::Int16
	lcf::Int16
	hcf::Int16
	lcs::Int16
	hcs::Int16
	year::Int16
	day::Int16
	hour::Int16
	minute::Int16
	sec::Int16
	timbas::Int16
	trwf::Int16
	grnors::Int16
	grnofr::Int16
	grnlof::Int16
	gaps::Int16
	otrav::Int16
	d1::Float32
	f1::Float32
	d2::Float32
	f2::Float32
	ungpow::Float32
	unscale::Float32
	ntr::Int32
	mark::Int16
	unass::Int16
end

segy_count = Dict{AbstractString,Int32}()
segy_count["tracl"]  = 0
segy_count["tracr"]  = 4
segy_count["fldr"]   = 8
segy_count["tracf"]  = 12
segy_count["ep"]     = 16
segy_count["cdp"]    = 20
segy_count["cdpt"]   = 24
segy_count["trid"]   = 28
segy_count["nva"]    = 30
segy_count["nhs"]    = 32
segy_count["duse"]   = 34
segy_count["offset"] = 36
segy_count["gelev"]  = 40
segy_count["selev"]  = 44
segy_count["sdepth"] = 48
segy_count["gdel"]   = 52
segy_count["sdel"]   = 56
segy_count["swdep"]  = 60
segy_count["gwdep"]  = 64
segy_count["scalel"] = 68
segy_count["scalco"] = 70
segy_count["sx"]     = 72
segy_count["sy"]     = 76
segy_count["gx"]     = 80
segy_count["gy"]     = 84
segy_count["counit"] = 88
segy_count["wevel"]  = 90
segy_count["swevel"] = 92
segy_count["sut"]    = 94
segy_count["gut"]    = 96
segy_count["sstat"]  = 98
segy_count["gstat"]  = 100
segy_count["tstat"]  = 102
segy_count["laga"]   = 104
segy_count["lagb"]   = 106
segy_count["delrt"]  = 108
segy_count["muts"]   = 110
segy_count["mute"]   = 112
segy_count["ns"]     = 114
segy_count["dt"]     = 116
segy_count["gain"]   = 118
segy_count["igc"]    = 120
segy_count["igi"]    = 122
segy_count["corr"]   = 124
segy_count["sfs"]    = 126
segy_count["sfe"]    = 128
segy_count["slen"]   = 130
segy_count["styp"]   = 132
segy_count["stas"]   = 134
segy_count["stae"]   = 136
segy_count["tatyp"]  = 138
segy_count["afilf"]  = 140
segy_count["afils"]  = 142
segy_count["nofilf"] = 144
segy_count["nofils"] = 146
segy_count["lcf"]    = 148
segy_count["hcf"]    = 150
segy_count["lcs"]    = 152
segy_count["hcs"]    = 154
segy_count["year"]   = 156
segy_count["day"]    = 158
segy_count["hour"]   = 160
segy_count["minute"] = 162
segy_count["sec"]    = 164
segy_count["timbas"] = 166
segy_count["trwf"]   = 168
segy_count["grnors"] = 170
segy_count["grnofr"] = 172
segy_count["grnlof"] = 174
segy_count["gaps"]   = 176
segy_count["otrav"]  = 178
segy_count["d1"]     = 180
segy_count["f1"]     = 184
segy_count["d2"]     = 188
segy_count["f2"]     = 192
segy_count["ungpow"] = 196
segy_count["unscale"]= 200
segy_count["ntr"]    = 204
segy_count["mark"]   = 208
segy_count["unass"]  = 210
segy_count["trace"]  = 240

function InitSegyHeader()
	h = SegyHeader(0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0)
	return h
end

function GrabSegyHeader(stream,swap_bytes,nt,file_header_size,j)
	position = file_header_size + (240+nt*4)*(j-1)
	seek(stream,position)
	h = InitSegyHeader()
	if (swap_bytes == false)
		h.tracl  = read(stream,typeof(h.tracl))
		h.tracr  = read(stream,typeof(h.tracr))
		h.fldr   = read(stream,typeof(h.fldr))
		h.tracf  = read(stream,typeof(h.tracf))
		h.ep     = read(stream,typeof(h.ep))
		h.cdp    = read(stream,typeof(h.cdp))
		h.cdpt   = read(stream,typeof(h.cdpt))
		h.trid   = read(stream,typeof(h.trid))
		h.nva    = read(stream,typeof(h.nva))
		h.nhs    = read(stream,typeof(h.nhs))
		h.duse   = read(stream,typeof(h.duse))
		h.offset = read(stream,typeof(h.offset))
		h.gelev  = read(stream,typeof(h.gelev))
		h.selev  = read(stream,typeof(h.selev))
		h.sdepth = read(stream,typeof(h.sdepth))
		h.gdel   = read(stream,typeof(h.gdel))
		h.sdel   = read(stream,typeof(h.sdel))
		h.swdep  = read(stream,typeof(h.swdep))
		h.gwdep  = read(stream,typeof(h.gwdep))
		h.scalel = read(stream,typeof(h.scalel))
		h.scalco = read(stream,typeof(h.scalco))
		h.sx     = read(stream,typeof(h.sx))
		h.sy     = read(stream,typeof(h.sy))
		h.gx     = read(stream,typeof(h.gx))
		h.gy     = read(stream,typeof(h.gy))
		h.counit = read(stream,typeof(h.counit))
		h.wevel  = read(stream,typeof(h.wevel))
		h.swevel = read(stream,typeof(h.swevel))
		h.sut    = read(stream,typeof(h.sut))
		h.gut    = read(stream,typeof(h.gut))
		h.sstat  = read(stream,typeof(h.sstat))
		h.gstat  = read(stream,typeof(h.gstat))
		h.tstat  = read(stream,typeof(h.tstat))
		h.laga   = read(stream,typeof(h.laga))
		h.lagb   = read(stream,typeof(h.lagb))
		h.delrt  = read(stream,typeof(h.delrt))
		h.muts   = read(stream,typeof(h.muts))
		h.mute   = read(stream,typeof(h.mute))
		h.ns     = read(stream,typeof(h.ns))
		h.dt     = read(stream,typeof(h.dt))
		h.gain   = read(stream,typeof(h.gain))
		h.igc    = read(stream,typeof(h.igc))
		h.igi    = read(stream,typeof(h.igi))
		h.corr   = read(stream,typeof(h.corr))
		h.sfs    = read(stream,typeof(h.sfs))
		h.sfe    = read(stream,typeof(h.sfe))
		h.slen   = read(stream,typeof(h.slen))
		h.styp   = read(stream,typeof(h.styp))
		h.stas   = read(stream,typeof(h.stas))
		h.stae   = read(stream,typeof(h.stae))
		h.tatyp  = read(stream,typeof(h.tatyp))
		h.afilf  = read(stream,typeof(h.afilf))
		h.afils  = read(stream,typeof(h.afils))
		h.nofilf = read(stream,typeof(h.nofilf))
		h.nofils = read(stream,typeof(h.nofils))
		h.lcf    = read(stream,typeof(h.lcf))
		h.hcf    = read(stream,typeof(h.hcf))
		h.lcs    = read(stream,typeof(h.lcs))
		h.hcs    = read(stream,typeof(h.hcs))
		h.year   = read(stream,typeof(h.year))
		h.day    = read(stream,typeof(h.day))
		h.hour   = read(stream,typeof(h.hour))
		h.minute = read(stream,typeof(h.minute))
		h.sec    = read(stream,typeof(h.sec))
		h.timbas = read(stream,typeof(h.timbas))
		h.trwf   = read(stream,typeof(h.trwf))
		h.grnors = read(stream,typeof(h.grnors))
		h.grnofr = read(stream,typeof(h.grnofr))
		h.grnlof = read(stream,typeof(h.grnlof))
		h.gaps   = read(stream,typeof(h.gaps))
		h.otrav  = read(stream,typeof(h.otrav))
		h.d1     = read(stream,typeof(h.d1))
		h.f1     = read(stream,typeof(h.f1))
		h.d2     = read(stream,typeof(h.d2))
		h.f2     = read(stream,typeof(h.f2))
		h.ungpow = read(stream,typeof(h.ungpow))
		h.unscale= read(stream,typeof(h.unscale))
		h.ntr    = read(stream,typeof(h.ntr))
		h.mark   = read(stream,typeof(h.mark))
		h.unass  = read(stream,typeof(h.unass))
	else
		h.tracl  = bswap(read(stream,typeof(h.tracl)))
		h.tracr  = bswap(read(stream,typeof(h.tracr)))
		h.fldr   = bswap(read(stream,typeof(h.fldr)))
		h.tracf  = bswap(read(stream,typeof(h.tracf)))
		h.ep     = bswap(read(stream,typeof(h.ep)))
		h.cdp    = bswap(read(stream,typeof(h.cdp)))
		h.cdpt   = bswap(read(stream,typeof(h.cdpt)))
		h.trid   = bswap(read(stream,typeof(h.trid)))
		h.nva    = bswap(read(stream,typeof(h.nva)))
		h.nhs    = bswap(read(stream,typeof(h.nhs)))
		h.duse   = bswap(read(stream,typeof(h.duse)))
		h.offset = bswap(read(stream,typeof(h.offset)))
		h.gelev  = bswap(read(stream,typeof(h.gelev)))
		h.selev  = bswap(read(stream,typeof(h.selev)))
		h.sdepth = bswap(read(stream,typeof(h.sdepth)))
		h.gdel   = bswap(read(stream,typeof(h.gdel)))
		h.sdel   = bswap(read(stream,typeof(h.sdel)))
		h.swdep  = bswap(read(stream,typeof(h.swdep)))
		h.gwdep  = bswap(read(stream,typeof(h.gwdep)))
		h.scalel = bswap(read(stream,typeof(h.scalel)))
		h.scalco = bswap(read(stream,typeof(h.scalco)))
		h.sx     = bswap(read(stream,typeof(h.sx)))
		h.sy     = bswap(read(stream,typeof(h.sy)))
		h.gx     = bswap(read(stream,typeof(h.gx)))
		h.gy     = bswap(read(stream,typeof(h.gy)))
		h.counit = bswap(read(stream,typeof(h.counit)))
		h.wevel  = bswap(read(stream,typeof(h.wevel)))
		h.swevel = bswap(read(stream,typeof(h.swevel)))
		h.sut    = bswap(read(stream,typeof(h.sut)))
		h.gut    = bswap(read(stream,typeof(h.gut)))
		h.sstat  = bswap(read(stream,typeof(h.sstat)))
		h.gstat  = bswap(read(stream,typeof(h.gstat)))
		h.tstat  = bswap(read(stream,typeof(h.tstat)))
		h.laga   = bswap(read(stream,typeof(h.laga)))
		h.lagb   = bswap(read(stream,typeof(h.lagb)))
		h.delrt  = bswap(read(stream,typeof(h.delrt)))
		h.muts   = bswap(read(stream,typeof(h.muts)))
		h.mute   = bswap(read(stream,typeof(h.mute)))
		h.ns     = bswap(read(stream,typeof(h.ns)))
		h.dt     = bswap(read(stream,typeof(h.dt)))
		h.gain   = bswap(read(stream,typeof(h.gain)))
		h.igc    = bswap(read(stream,typeof(h.igc)))
		h.igi    = bswap(read(stream,typeof(h.igi)))
		h.corr   = bswap(read(stream,typeof(h.corr)))
		h.sfs    = bswap(read(stream,typeof(h.sfs)))
		h.sfe    = bswap(read(stream,typeof(h.sfe)))
		h.slen   = bswap(read(stream,typeof(h.slen)))
		h.styp   = bswap(read(stream,typeof(h.styp)))
		h.stas   = bswap(read(stream,typeof(h.stas)))
		h.stae   = bswap(read(stream,typeof(h.stae)))
		h.tatyp  = bswap(read(stream,typeof(h.tatyp)))
		h.afilf  = bswap(read(stream,typeof(h.afilf)))
		h.afils  = bswap(read(stream,typeof(h.afils)))
		h.nofilf = bswap(read(stream,typeof(h.nofilf)))
		h.nofils = bswap(read(stream,typeof(h.nofils)))
		h.lcf    = bswap(read(stream,typeof(h.lcf)))
		h.hcf    = bswap(read(stream,typeof(h.hcf)))
		h.lcs    = bswap(read(stream,typeof(h.lcs)))
		h.hcs    = bswap(read(stream,typeof(h.hcs)))
		h.year   = bswap(read(stream,typeof(h.year)))
		h.day    = bswap(read(stream,typeof(h.day)))
		h.hour   = bswap(read(stream,typeof(h.hour)))
		h.minute = bswap(read(stream,typeof(h.minute)))
		h.sec    = bswap(read(stream,typeof(h.sec)))
		h.timbas = bswap(read(stream,typeof(h.timbas)))
		h.trwf   = bswap(read(stream,typeof(h.trwf)))
		h.grnors = bswap(read(stream,typeof(h.grnors)))
		h.grnofr = bswap(read(stream,typeof(h.grnofr)))
		h.grnlof = bswap(read(stream,typeof(h.grnlof)))
		h.gaps   = bswap(read(stream,typeof(h.gaps)))
		h.otrav  = bswap(read(stream,typeof(h.otrav)))
		h.d1     = bswap(read(stream,typeof(h.d1)))
		h.f1     = bswap(read(stream,typeof(h.f1)))
		h.d2     = bswap(read(stream,typeof(h.d2)))
		h.f2     = bswap(read(stream,typeof(h.f2)))
		h.ungpow = bswap(read(stream,typeof(h.ungpow)))
		h.unscale= bswap(read(stream,typeof(h.unscale)))
		h.ntr    = bswap(read(stream,typeof(h.ntr)))
		h.mark   = bswap(read(stream,typeof(h.mark)))
		h.unass  = bswap(read(stream,typeof(h.unass)))
	end

	return h
end

function PutSegyHeader(stream,h,nt,file_header_size,j)
	position = file_header_size + (240+nt*4)*(j-1)
	seek(stream,position)
	write(stream,h.tracl)
	write(stream,h.tracr)
	write(stream,h.fldr)
	write(stream,h.tracf)
	write(stream,h.ep)
	write(stream,h.cdp)
	write(stream,h.cdpt)
	write(stream,h.trid)
	write(stream,h.nva)
	write(stream,h.nhs)
	write(stream,h.duse)
	write(stream,h.offset)
	write(stream,h.gelev)
	write(stream,h.selev)
	write(stream,h.sdepth)
	write(stream,h.gdel)
	write(stream,h.sdel)
	write(stream,h.swdep)
	write(stream,h.gwdep)
	write(stream,h.scalel)
	write(stream,h.scalco)
	write(stream,h.sx)
	write(stream,h.sy)
	write(stream,h.gx)
	write(stream,h.gy)
	write(stream,h.counit)
	write(stream,h.wevel)
	write(stream,h.swevel)
	write(stream,h.sut)
	write(stream,h.gut)
	write(stream,h.sstat)
	write(stream,h.gstat)
	write(stream,h.tstat)
	write(stream,h.laga)
	write(stream,h.lagb)
	write(stream,h.delrt)
	write(stream,h.muts)
	write(stream,h.mute)
	write(stream,h.ns)
	write(stream,h.dt)
	write(stream,h.gain)
	write(stream,h.igc)
	write(stream,h.igi)
	write(stream,h.corr)
	write(stream,h.sfs)
	write(stream,h.sfe)
	write(stream,h.slen)
	write(stream,h.styp)
	write(stream,h.stas)
	write(stream,h.stae)
	write(stream,h.tatyp)
	write(stream,h.afilf)
	write(stream,h.afils)
	write(stream,h.nofilf)
	write(stream,h.nofils)
	write(stream,h.lcf)
	write(stream,h.hcf)
	write(stream,h.lcs)
	write(stream,h.hcs)
	write(stream,h.year)
	write(stream,h.day)
	write(stream,h.hour)
	write(stream,h.minute)
	write(stream,h.sec)
	write(stream,h.timbas)
	write(stream,h.trwf)
	write(stream,h.grnors)
	write(stream,h.grnofr)
	write(stream,h.grnlof)
	write(stream,h.gaps)
	write(stream,h.otrav)
	write(stream,h.d1)
	write(stream,h.f1)
	write(stream,h.d2)
	write(stream,h.f2)
	write(stream,h.ungpow)
	write(stream,h.unscale)
	write(stream,h.ntr)
	write(stream,h.mark)
	write(stream,h.unass)
end

function MapHeaders(h_in,j,map_type)

	if map_type=="SegyToSeis"
		#scalco = abs(convert(Float32,h_in[1].scalco)) < 10 ? convert(Float32,h_in[1].scalco) : sign(convert(Float32,h_in[1].scalco))*log10(abs(convert(Float32,h_in[1].scalco)))
		scalco = sign(convert(Float32,h_in[1].scalco))*log10(abs(convert(Float32,h_in[1].scalco)))

		#scalel = abs(convert(Float32,h_in[1].scalel)) < 10 ? convert(Float32,h_in[1].scalel) : sign(convert(Float32,h_in[1].scalco))*log10(abs(convert(Float32,h_in[1].scalco)))
		scalel = sign(convert(Float32,h_in[1].scalco))*log10(abs(convert(Float32,h_in[1].scalco)))
		h_out = InitSeisHeader()
		h_out.tracenum = convert(typeof(h_out.tracenum),j)
		h_out.o1 = convert(typeof(h_out.o1),0)
		h_out.n1 = convert(typeof(h_out.n1),h_in[1].ns)
		h_out.d1 = convert(typeof(h_out.d1),h_in[1].dt/1000000)
		h_out.sx = scalco >= 0 ? convert(typeof(h_out.sx),h_in[1].sx)*10^scalco : convert(typeof(h_out.sx),h_in[1].sx)/(10^abs(scalco))
		h_out.sy = scalco >= 0 ? convert(typeof(h_out.sy),h_in[1].sy)*10^scalco : convert(typeof(h_out.sx),h_in[1].sy)/(10^abs(scalco))
		h_out.gx = scalco >= 0 ? convert(typeof(h_out.gx),h_in[1].gx)*10^scalco : convert(typeof(h_out.gx),h_in[1].gx)/(10^abs(scalco))
		h_out.gy = scalco >= 0 ? convert(typeof(h_out.gy),h_in[1].gy)*10^scalco : convert(typeof(h_out.gy),h_in[1].gy)/(10^abs(scalco))
		h_out.h = convert(typeof(h_out.h),h_in[1].offset)
		h_out.isx = convert(typeof(h_out.isx),h_in[1].ep)
		h_out.imx = convert(typeof(h_out.imx),h_in[1].cdp)
		h_out.selev = scalel >= 0 ? convert(typeof(h_out.selev),h_in[1].selev)*(10^scalel) : convert(typeof(h_out.selev),h_in[1].selev)/(10^abs(scalel))
		h_out.gelev = scalel >= 0 ? convert(typeof(h_out.gelev),h_in[1].gelev)*(10^scalel) : convert(typeof(h_out.gelev),h_in[1].gelev)/(10^abs(scalel))
		h_out.trid = convert(typeof(h_out.trid),h_in[1].trid)
	else
		h_out = InitSegyHeader()
		h_out.tracl = convert(typeof(h_out.tracl),j)
		h_out.tracr = convert(typeof(h_out.tracr),j)
		h_out.scalel = convert(typeof(h_out.scalel),-3)
		h_out.scalco = convert(typeof(h_out.scalco),-3)
		h_out.counit = convert(typeof(h_out.counit),1)
		h_out.gain = convert(typeof(h_out.gain),1)
		h_out.corr = convert(typeof(h_out.corr),1)
		h_out.ns = convert(typeof(h_out.ns),Int(h_in[1].n1))
		h_out.dt = convert(typeof(h_out.dt),round(Int,h_in[1].d1*1000000))
		h_out.sx = convert(typeof(h_out.sx),round(Int,h_in[1].sx*1000))
		h_out.sy = convert(typeof(h_out.sy),round(Int,h_in[1].sy*1000))
		h_out.gx = convert(typeof(h_out.gx),round(Int,h_in[1].gx*1000))
		h_out.gy = convert(typeof(h_out.gy),round(Int,h_in[1].gy*1000))
		h_out.offset = convert(typeof(h_out.offset),round(Int,h_in[1].h))
		h_out.ep = convert(typeof(h_out.ep),Int(h_in[1].isx))
		h_out.cdp = convert(typeof(h_out.cdp),Int(h_in[1].imx))
		h_out.selev = convert(typeof(h_out.selev),round(Int,h_in[1].selev*1000))
		h_out.gelev = convert(typeof(h_out.gelev),round(Int,h_in[1].gelev*1000))
		h_out.trid = convert(typeof(h_out.trid),h_in[1].trid)
	end

	return h_out
end

import Base.convert

bitstype 32 IBMFloat32

ieeeOfPieces(fr::UInt32, exp::Int32, sgn::UInt32) = reinterpret(Float32, convert(UInt32,fr >>> 9) | convert(UInt32,exp << 23) | sgn) :: Float32
import Base.convert

function convert(::Type{Float32}, ibm::IBMFloat32)
  local fr::UInt32 = ntoh(reinterpret(UInt32, ibm))
  local sgn::UInt32 = fr & 0x80000000 # save sign
  fr <<= 1 # shift sign out
  local exp::Int32 = convert(Int32,fr >>> 25) # save exponent
  fr <<= 7 # shift exponent out

  if (fr == convert(UInt32,0))
    zero(Float32)
  else
    # normalize the signficand
    local norm::UInt32 = leading_zeros(fr)
    fr <<= norm
    exp = (exp << 2) - 130 - norm

    # exp <= 0 --> ieee(0,0,sgn)
    # exp >= 255 --> ieee(0,255,sgn)
    # else -> ieee(fr<<1, exp, sgn)
    local clexp::Int32 = exp & convert(Int32,0xFF)
    ieeeOfPieces(clexp == exp ? fr << 1 : convert(UInt32,0), clexp, sgn)
  end
end

