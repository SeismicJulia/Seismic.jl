type SepShot3C
	ux
	uy
	uz
	up
	us
	vp
	vs
	fmin
	fmax
	decomp
	H
	verbose
end

function Sep1Shot(shot)

	run( `wavesep H=${shot.H}
	ux=${shot.ux} uy=${shot.uy} uz=${shot.uz} 
	up=${shot.up} us=${shot.us} 
	vp=${shot.vp} vs=${shot.vs}
	fmin=${shot.fmin} fmax=${shot.fmax}
	verbose=${shot.verbose}`)
	return(1)

end

function WaveSep(ux,uy,uz,up,us,param)

	H = get(param,"H","y") # flag for exact Helmholtz decomp op (y), or adjoint of inverse Helmholtz decomp op (n)
	decomp = get(param,"decomp","y") # flag for decomposition from data components to potentials (y), or recomposition from potentials to data components (n)
	vp = get(param,"vp",2000) # the P-wave velocity 
	vs = get(param,"vs",1000) # the S-wave velocity 
	fmin = get(param,"fmin",0)  # min frequency to process (Hz)
	fmax = get(param,"fmax",80) # max frequency to process (Hz)
	verbose = get(param,"verbose","n") # flag for error / debugging messages
	isx = get(param,"isx",[0]) # array of source bin numbers in the x direction
	isy = get(param,"isy",[0]) # array of source bin numbers in the y direction
	nshot = length(isx)	

	shot_list = Array(SepShot3C,nshot)
	for ishot = 1 : nshot
		shot_list[ishot] = SepShot3C(0,0,0,0,0,0,0,0,0,0,0,0)
		shot_list[ishot].ux = join([ux "_shot_" int(floor(isx[ishot])) "_" int(floor(isy[ishot]))])
		shot_list[ishot].uy = join([uy "_shot_" int(floor(isx[ishot])) "_" int(floor(isy[ishot]))])
		shot_list[ishot].uz = join([uz "_shot_" int(floor(isx[ishot])) "_" int(floor(isy[ishot]))])
		shot_list[ishot].up = join([up "_shot_" int(floor(isx[ishot])) "_" int(floor(isy[ishot]))])
		shot_list[ishot].us = join([us "_shot_" int(floor(isx[ishot])) "_" int(floor(isy[ishot]))])
		shot_list[ishot].vp = vp
		shot_list[ishot].vs = vs
		shot_list[ishot].fmin = fmin
		shot_list[ishot].fmax = fmax
		shot_list[ishot].decomp = decomp
		shot_list[ishot].H = H
		shot_list[ishot].verbose = verbose
	end
	if (decomp == "y")
		for ishot = 1 : nshot
			SeisWindow(ux,shot_list[ishot].ux,["isx","isy"],[isx[ishot],isy[ishot]],[isx[ishot],isy[ishot]])
			SeisWindow(uy,shot_list[ishot].uy,["isx","isy"],[isx[ishot],isy[ishot]],[isx[ishot],isy[ishot]])
			SeisWindow(uz,shot_list[ishot].uz,["isx","isy"],[isx[ishot],isy[ishot]],[isx[ishot],isy[ishot]])
		end
		a = pmap(Sep1Shot,shot_list)
		j = 1    
		for ishot = 1 : nshot
			up_shot,h_shot = SeisRead(shot_list[ishot].up)
			for itrace = 1 : size(up_shot,2)
				h_shot[ishot].tracenum = j + itrace	
			end
			SeisWrite(up,up_shot,h_shot,["itrace"=>j])
			us_shot,h_shot = SeisRead(shot_list[ishot].us)
			for itrace = 1 : size(us_shot,2)
				h_shot[ishot].tracenum = j + itrace	
			end
			SeisWrite(us,us_shot,h_shot,["itrace"=>j])
			j += size(up_shot,2)
			SeisRemove(shot_list[ishot].up)
			SeisRemove(shot_list[ishot].us)
			SeisRemove(shot_list[ishot].ux)
			SeisRemove(shot_list[ishot].uy)
			SeisRemove(shot_list[ishot].uz)

		end	
	else
		for ishot = 1 : nshot
			SeisWindow(up,shot_list[ishot].up,["isx","isy"],[isx[ishot],isy[ishot]],[isx[ishot],isy[ishot]])
			SeisWindow(us,shot_list[ishot].us,["isx","isy"],[isx[ishot],isy[ishot]],[isx[ishot],isy[ishot]])
		end
		a = pmap(Sep1Shot,shot_list)
		j = 1    
		for ishot = 1 : nshot
			ux_shot,h_shot = SeisRead(shot_list[ishot].ux)
			SeisWrite(ux,ux_shot,h_shot,["itrace"=>j])
			uy_shot,h_shot = SeisRead(shot_list[ishot].uy)
			SeisWrite(uy,uy_shot,h_shot,["itrace"=>j])
			uz_shot,h_shot = SeisRead(shot_list[ishot].uz)
			SeisWrite(uz,uz_shot,h_shot,["itrace"=>j])
			j += size(ux_shot,2)
			SeisRemove(shot_list[ishot].up)
			SeisRemove(shot_list[ishot].us)
			SeisRemove(shot_list[ishot].ux)
			SeisRemove(shot_list[ishot].uy)
			SeisRemove(shot_list[ishot].uz)
		end	
	end
end
