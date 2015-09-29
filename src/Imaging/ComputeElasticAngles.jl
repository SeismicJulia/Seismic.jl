function ComputeElasticAngles(angpx,angpy,angsx,angsy,param)

	vp = get(param,"vp","vp") # seis file containing the p-wave velocity (should have same x, y, and z dimensions as the desired image)
	vs = get(param,"vs","vs") # seis file containing the s-wave velocity (should have same x, y, and z dimensions as the desired image)
	dip_flag = get(param,"dip_flag","n") # flag for correcting angles to be relative to reflector normal
	dipx_name = get(param,"dipx","NULL") # seis file containing the reflector normals relative to vertical in the x direction (should have same x, y, and z dimensions as the desired image)
	dipy_name = get(param,"dipy","NULL") # seis file containing the reflector normals relative to vertical in the y direction (should have same x, y, and z dimensions as the desired image)
	wav = get(param,"wav","wav") # seis file containing the source wavelet (in time domain) 
	sz = get(param,"sz",0) # source depth (Dev: read this from source wavelet file for variable source depth)
	gz = get(param,"gz",0) # receiver depth (Dev: read this from data file for variable source depth (but then what to do in fwd op?))
	nhx = get(param,"nhx",101)  # number of offset bins
	ohx = get(param,"ohx",1000) # min offset (surface offset in the data)
	dhx = get(param,"dhx",10)   # offset increment
	nhy = get(param,"nhy",101)  # number of offset bins
	ohy = get(param,"ohy",1000) # min offset (surface offset in the data)
	dhy = get(param,"dhy",10)   # offset increment
	pade_flag = get(param,"pade_flag","n") # flag for Pade Fourier correction
	fmin = get(param,"fmin",0)  # min frequency to process (Hz)
	fmax = get(param,"fmax",80) # max frequency to process (Hz)
	padt = get(param,"padt",2)  # pad factor for the time axis
	padx = get(param,"padx",2) # pad factor for the spatial axes
	omp = get(param,"omp",1) # number of shared memory threads to use for frequency slice processing 
	verbose = get(param,"verbose","n") # flag for error / debugging messages
	sx = get(param,"sx",[0]) # array of source X positions (meters)
	sy = get(param,"sy",[0]) # array of source Y positions (meters)
	nshot = length(sx)	

	if (dip_flag=="y")
		dipx,h = SeisRead(dipx_name)
		dipy,h = SeisRead(dipy_name)
	end	

	v,h = SeisRead(vp)
	min_imx = h[1].imx
	max_imx = h[end].imx
	dmx = dhx
	nmx = max_imx - min_imx + 1
	min_imy = h[1].imy
	max_imy = h[end].imy
	dmy = dhy
	nmy = max_imy - min_imy + 1
	nz = h[1].n1
	dz = h[1].d1
	oz = h[1].o1
	min_gx = h[1].mx
	max_gx = h[end].mx
	min_gy = h[1].my
	max_gy = h[end].my

	shot_list = Array(ElasticShotAngles,nshot)
	for ishot = 1 : nshot
		shot_list[ishot] = ElasticShotAngles(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		shot_list[ishot].angpx = join(["angpx_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].angpy = join(["angpy_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].angsx = join(["angsx_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].angsy = join(["angsy_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].vp = vp
		shot_list[ishot].vs = vs
		shot_list[ishot].dipx = dipx_name
		shot_list[ishot].dipy = dipy_name
		shot_list[ishot].wav = wav
		shot_list[ishot].sx = sx[ishot]
		shot_list[ishot].sy = sy[ishot]
		shot_list[ishot].sz = sz
		shot_list[ishot].ohx = ohx
		shot_list[ishot].dhx = dhx
		shot_list[ishot].nhx = nhx
		shot_list[ishot].ohy = ohy
		shot_list[ishot].dhy = dhy
		shot_list[ishot].nhy = nhy
		shot_list[ishot].fmin = fmin
		shot_list[ishot].fmax = fmax
		shot_list[ishot].omp = omp
		shot_list[ishot].pade_flag = pade_flag
		shot_list[ishot].dip_flag = dip_flag
		shot_list[ishot].verbose = verbose	
	end

	a = pmap(compute_eangles,shot_list)

	j = 1
	for ishot = 1 : nshot
		angpx_shot,h_shot = SeisRead(shot_list[ishot].angpx)
		angpy_shot,h_shot = SeisRead(shot_list[ishot].angpy)
		angsx_shot,h_shot = SeisRead(shot_list[ishot].angsx)
		angsy_shot,h_shot = SeisRead(shot_list[ishot].angsy)
		SeisWrite(angpx,angpx_shot,h_shot,["itrace"=>j])
		SeisWrite(angpy,angpy_shot,h_shot,["itrace"=>j])
		SeisWrite(angsx,angsx_shot,h_shot,["itrace"=>j])
		SeisWrite(angsy,angsy_shot,h_shot,["itrace"=>j])
		j += size(angpx_shot,2)
		SeisRemove(shot_list[ishot].angpx)
		SeisRemove(shot_list[ishot].angpy)
		SeisRemove(shot_list[ishot].angsx)
		SeisRemove(shot_list[ishot].angsy)
	end

end

type ElasticShotAngles
	angpx
	angpy
	angsx
	angsy
	vp
	vs
	dipx
	dipy
	wav
	sx
	sy
	sz
	ohx
	dhx
	nhx
	ohy
	dhy
	nhy
	fmin
	fmax
	pade_flag
	dip_flag
	omp
	verbose
end

function compute_eangles(shot)

	ENV["OMP_NUM_THREADS"] = shot.omp
	if (find_library(["compute_eangles"]) == "")
		error("Couldn't find shared library compute_eangles.so, make sure it is installed and check LD_LIBRARY_PATH environmental variable.")
	end
	argv = ["0",
	join(["angpx=",shot.angpx]), join(["angpy=",shot.angpy]), join(["angsx=",shot.angsx]), join(["angsy=",shot.angsy]),
	join(["vp=",shot.vp]), join(["vs=",shot.vs]), 
	join(["dipx=",shot.dipx]), join(["dipy=",shot.dipy]), join(["wav=",shot.wav]),  
	join(["sx=",shot.sx]), join(["sy=",shot.sy]),  join(["sz=",shot.sz]),
	join(["ohx=",shot.ohx]),  join(["dhx=",shot.dhx]),  join(["nhx=",shot.nhx]), 
	join(["ohy=",shot.ohy]),  join(["dhy=",shot.dhy]),  join(["nhy=",shot.nhy]),   
	join(["fmin=",shot.fmin]),  join(["fmax=",shot.fmax]), 
	join(["pade_flag=",shot.pade_flag]),  join(["dip_flag=",shot.dip_flag]),  join(["verbose=",shot.verbose]) ] 
	a = ccall((:main, "compute_eangles"), Int32, (Int32, Ptr{Ptr{Uint8}}), length(argv), argv)
	return(a)


end

