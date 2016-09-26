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
    
    run( `wavesep H=$(shot.H)
        ux=$(shot.ux) uy=$(shot.uy) uz=$(shot.uz) 
	up=$(shot.up) us=$(shot.us) 
	vp=$(shot.vp) vs=$(shot.vs)
	fmin=$(shot.fmin) fmax=$(shot.fmax)
	verbose=$(shot.verbose)`)
    return(1)
    
end

function WaveSep(ux, uy, uz, up, us;
                 H=true, decomp=true, vp=2000, vs=1000, fmin=0, fmax=80,
                 verbose=true, isx=0, isy=0)

    # H :  flag for exact Helmholtz decomp op (true), or adjoint of inverse
    # Helmholtz decomp op (false)
    # decomp : flag for decomposition from data components to potentials (true),
    # or recomposition from potentials to data components (false)
    # vp :  the P-wave velocity 
    # vs : the S-wave velocity 
    # fmin : min frequency to process (Hz)
    # fmax : max frequency to process (Hz)
    # verbose : flag for error / debugging messages
    # isx : array of source bin numbers in the x direction
    # isy : array of source bin numbers in the y direction
    
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
	shot_list[ishot].decomp = decomp == true ? "y" : "n"
	shot_list[ishot].H = H == true ? "y" : "n"
	shot_list[ishot].verbose = verbose == true ? "y" : "n"
    end
    if (decomp == "y")
	for ishot = 1 : nshot
	    SeisWindow(ux, shot_list[ishot].ux, key=["isx","isy"],
                       minval=[isx[ishot],isy[ishot]],
                       maxval=[isx[ishot],isy[ishot]])
	    SeisWindow(uy, shot_list[ishot].uy, key=["isx","isy"],
                       minval=[isx[ishot],isy[ishot]],
                       maxval=[isx[ishot],isy[ishot]])
	    SeisWindow(uz, shot_list[ishot].uz, key=["isx","isy"],
                       minval=[isx[ishot],isy[ishot]],
                       maxval=[isx[ishot],isy[ishot]])
	end
	a = pmap(Sep1Shot,shot_list)
	j = 1    
	for ishot = 1 : nshot
	    up_shot,h_shot,e = SeisRead(shot_list[ishot].up)
	    for itrace = 1 : size(up_shot,2)
		h_shot[ishot].tracenum = j + itrace	
	    end
	    SeisWrite(up,up_shot,h_shot,itrace=j)
	    us_shot,h_shot,e = SeisRead(shot_list[ishot].us)
	    for itrace = 1 : size(us_shot,2)
		h_shot[ishot].tracenum = j + itrace	
	    end
	    SeisWrite(us,us_shot,h_shot,itrace=j)
	    j += size(up_shot,2)
	    SeisRemove(shot_list[ishot].up)
	    SeisRemove(shot_list[ishot].us)
	    SeisRemove(shot_list[ishot].ux)
	    SeisRemove(shot_list[ishot].uy)
	    SeisRemove(shot_list[ishot].uz)
            
	end	
    else
	for ishot = 1 : nshot
	    SeisWindow(up, shot_list[ishot].up, key=["isx","isy"],
                       minval=[isx[ishot],isy[ishot]],
                       maxval=[isx[ishot],isy[ishot]])
	    SeisWindow(us, shot_list[ishot].us, key=["isx","isy"],
                       minval=[isx[ishot],isy[ishot]],
                       maxval=[isx[ishot],isy[ishot]])
	end
	a = pmap(Sep1Shot,shot_list)
	j = 1    
	for ishot = 1 : nshot
	    ux_shot,h_shot,e = SeisRead(shot_list[ishot].ux)
	    SeisWrite(ux,ux_shot,h_shot,itrace=j)
	    uy_shot,h_shot,e = SeisRead(shot_list[ishot].uy)
	    SeisWrite(uy,uy_shot,h_shot,itrace=j)
	    uz_shot,h_shot,e = SeisRead(shot_list[ishot].uz)
	    SeisWrite(uz,uz_shot,h_shot,itrace=j)
	    j += size(ux_shot,2)
	    SeisRemove(shot_list[ishot].up)
	    SeisRemove(shot_list[ishot].us)
	    SeisRemove(shot_list[ishot].ux)
	    SeisRemove(shot_list[ishot].uy)
	    SeisRemove(shot_list[ishot].uz)
	end	
    end
end
