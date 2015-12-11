function ShotProfileWEM(m::ASCIIString,d::ASCIIString,param::Dict{Any,Any})
# 
# ShotProfileWEM: Shot Profile Wave Equation Migration and Demigration 
# of 3D isotropic data.
#
# IN   
#      adjoint operator (param["adj"] = true)
#              d : data file name
#      forward operator (param["adj"] = false)
#              m : image file name
#
#      param={"adj"=>true,
#			  "vel"=>"vel",
#             "angx"=>"angx","angy"=>"angy",
#             "wav"=>"wav",
#             "sz"=>sz,"gz"=>gz,
#             "nhx"=>nhx,"ohx"=>ohx,"dhx"=>dhx,"nhy"=>nhy,"ohy"=>ohy,"dhy"=>dhy,
#             "pade_flag"=>"n"
#
#      Where the parameters stored in param are:
#
#      adj : flag for adjoint (migration), or forward (demigration)
#      vel : seis file containing the velocity (should have same x, y and z dimensions as the desired image)
#      angx : seis file containing incidence angles in the x direction for each shot
#      angy : seis file containing incidence angles in the y direction for each shot
#      wav : seis file containing the source wavelet (in time domain) 
#      sz : source depth (Dev: read this from source wavelet file for variable source depth)
#      gz : receiver depth (Dev: read this from data file for variable source depth (but then what to do in fwd op?))
#      nhx : number of offset bins
#      ohx : min offset (surface offset of the data)
#      dhx : offset increment
#      nhy : number of offset bins
#      ohy : min offset (surface offset of the data)
#      dhy : offset increment
#      pade_flag : flag for Pade Fourier correction
#      nangx : number of angle bins in x direction
#      oangx : min angle in x direction (angle between source incidence angle and reflector normal in Degrees)
#      dangx : angle increment in x direction
#      nangy : number of angle bins in y direction
#      oangy : min angle in y direction (angle between source incidence angle and reflector normal in Degrees)
#      dangy : angle increment in y direction
#      fmin : min frequency to process (Hz)
#      fmax : max frequency to process (Hz)
#      padt : pad factor for the time axis
#      padx : pad factor for the spatial axes
#      verbose : flag for error / debugging messages
#      sx : array of source X positions (meters)
#      sy : array of source Y positions (meters)
# 
# OUT  note: outputs are specified in the function parameters with inputs.
#      if param["adj"] == true
#              mpp : PP image file name
#              mps : PS image file name
#              msp : SP image file name
#              mss : SS image file name
#      else
#              ux : data X-component file name
#              uy : data Y-component file name
#              uz : data Z-component file name
#
# Example of forward operator: param["adj"] = false
#
#	Make reflectivity and velocity models
# 	param = {"ot"=>0,"dt"=>10,"nt"=>500,"ox1"=>0,"dx1"=>10,"nx1"=>500,"tau1"=>a,"v1"=>a.*0 + 3,"amp"=>[1,-1,1,-1,1,-1,1,-1,1,-1,1,-1],"f0"=>a.*0 + 0.01,"sinusoidal"=>true,"sinusoidalA"=>200};
#	m = SeisLinearEvents(param); m = m[:,:];
#	v = m.*0 + 2000.;
#   
# Example with parameters provided by the user 
# 
#      param=["dt"=>0.004,"fc"=>20]
#      (t,w) = ricker(param)
#      plot(t,w) 
#      xlabel('Time (s)')

	adj = get(param,"adj",true) # flag for adjoint (migration), or forward (demigration)
	damping = get(param,"damping",1000.)  # damping for deconvolution imaging condition
	vel = get(param,"vel","vel") # seis file containing the velocity (should have same x and z dimensions as the desired image)
	angx = get(param,"angx","angx") # seis file containing incidence angles in the x direction for each shot
	angy = get(param,"angy","angy") # seis file containing incidence angles in the y direction for each shot
	wav = get(param,"wav","wav") # seis file containing the source wavelet (in time domain) 
	sz = get(param,"sz",0.) # source depth (Dev: read this from source wavelet file for variable source depth)
	gz = get(param,"gz",0.) # receiver depth (Dev: read this from data file for variable source depth (but then what to do in fwd op?))
	nhx = get(param,"nhx",101)  # number of offset bins
	ohx = get(param,"ohx",1000.) # min offset (surface offset in the data)
	dhx = get(param,"dhx",10.)   # offset increment
	nhy = get(param,"nhy",101)  # number of offset bins
	ohy = get(param,"ohy",1000.) # min offset (surface offset in the data)
	dhy = get(param,"dhy",10.)   # offset increment
	pade_flag = get(param,"pade_flag",false) # flag for Pade Fourier correction
	nangx = get(param,"nangx",1) # number of angle bins in x direction
	oangx = get(param,"oangx",0.) # min angle in x direction (angle between source incidence angle and reflector normal in Degrees)
	dangx = get(param,"dangx",1.) # angle increment in x direction
	nangy = get(param,"nangy",1) # number of angle bins in y direction
	oangy = get(param,"oangy",0.) # min angle in y direction (angle between source incidence angle and reflector normal in Degrees)
	dangy = get(param,"dangy",1.) # angle increment in y direction
	fmin = get(param,"fmin",0.)  # min frequency to process (Hz)
	fmax = get(param,"fmax",80.) # max frequency to process (Hz)
	padt = get(param,"padt",2)  # pad factor for the time axis
	padx = get(param,"padx",2) # pad factor for the spatial axes
	verbose = get(param,"verbose",false) # flag for error / debugging messages
	sx = get(param,"sx",[0.]) # array of source X positions (meters)
	sy = get(param,"sy",[0.]) # array of source Y positions (meters)
	nshot = length(sx)	

	v,h = SeisRead(vel)
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

	shot_list = Array(Shot,nshot)
	for ishot = 1 : nshot
		shot_list[ishot] = Shot(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		shot_list[ishot].m = join([m "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].d = join([d "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].angx = join([angx "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].angy = join([angy "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].vel = vel
		shot_list[ishot].wav = wav
		shot_list[ishot].damping = damping
		shot_list[ishot].sx = sx[ishot]
		shot_list[ishot].sy = sy[ishot]
		shot_list[ishot].sz = sz
		shot_list[ishot].gz = gz
		shot_list[ishot].ohx = ohx
		shot_list[ishot].dhx = dhx
		shot_list[ishot].nhx = nhx
		shot_list[ishot].ohy = ohy
		shot_list[ishot].dhy = dhy
		shot_list[ishot].nhy = nhy
		shot_list[ishot].fmin = fmin
		shot_list[ishot].fmax = fmax
		shot_list[ishot].padt = padt
		shot_list[ishot].padx = padx		
		shot_list[ishot].adj = (adj == true) ? "y" : "n"
		shot_list[ishot].pade_flag = (pade_flag == true) ? "y" : "n"
		shot_list[ishot].verbose = (verbose == true) ? "y" : "n"
	end

	if (nangx != 1 || nangy != 1)
		for ishot = 1 : nshot
			SeisWindow(angx,shot_list[ishot].angx,["key"=>["sx","sy"],"minval"=>[sx[ishot],sy[ishot]],"maxval"=>[sx[ishot],sy[ishot]]])
			SeisWindow(angy,shot_list[ishot].angy,["key"=>["sx","sy"],"minval"=>[sx[ishot],sy[ishot]],"maxval"=>[sx[ishot],sy[ishot]]])
		end
	end

	if (adj == true)
		for ishot = 1 : nshot
			SeisWindow(d,shot_list[ishot].d,["key"=>["sx","sy","gx","gy"],"minval"=>[sx[ishot],sy[ishot],min_gx,min_gy],"maxval"=>[sx[ishot],sy[ishot],max_gx,max_gy]])
		end
		a = pmap(shotwem,shot_list)
		j = 1    
		gather = zeros(Float32,nz,nangx*nangy)

		for imx = 1 : nmx
			for imy = 1 : nmy
				h = Header[]
				for iangx = 1 : nangx
					for iangy = 1 : nangy
						push!(h,Seismic.InitSeisHeader())		
						h[(iangx-1)*nangy + iangy].tracenum = convert(typeof(h[1].tracenum),j + (iangx-1)*nangy + iangy)
						h[(iangx-1)*nangy + iangy].o1 = convert(typeof(h[1].o1),0)
						h[(iangx-1)*nangy + iangy].n1 = convert(typeof(h[1].n1),nz)
						h[(iangx-1)*nangy + iangy].d1 = convert(typeof(h[1].d1),dz)
						h[(iangx-1)*nangy + iangy].imx = convert(typeof(h[1].imx),imx-1 + min_imx)
						h[(iangx-1)*nangy + iangy].imy = convert(typeof(h[1].imy),imy-1 + min_imy)			
						h[(iangx-1)*nangy + iangy].mx = convert(typeof(h[1].mx),dhx*(imx-1 + min_imx))
						h[(iangx-1)*nangy + iangy].my = convert(typeof(h[1].my),dhy*(imy-1 + min_imy))
						h[(iangx-1)*nangy + iangy].iang = convert(typeof(h[1].iang),iangx-1)		
						h[(iangx-1)*nangy + iangy].ang = convert(typeof(h[1].ang),(iangx-1)*dangx + oangx)
						h[(iangx-1)*nangy + iangy].iaz = convert(typeof(h[1].iaz),iangy-1)
						h[(iangx-1)*nangy + iangy].az = convert(typeof(h[1].az),(iangy-1)*dangy + oangy)
					end
				end
				SeisWrite(m,gather,h,["itrace"=>j])
				j += nangx*nangy
			end
		end

		m_m = join([m ".seisd"])
		m_h = join([m ".seish"])
		stream_m = open(m_m,"a+")
		stream_h = open(m_h,"a+")
		for ishot = 1 : nshot

			# println("reading ",shot_list[ishot].m)
			m_shot,h_shot  = SeisRead(shot_list[ishot].m)
			if (nangx != 1 || nangy != 1)
				angx_shot,h_ang = SeisRead(shot_list[ishot].angx)
				angy_shot,h_ang = SeisRead(shot_list[ishot].angy)
			else
				angx_shot = 0.*m_shot
				angy_shot = 0.*m_shot
			end	
			nx_shot = size(m_shot,2)
			for ix = 1 : nx_shot
				itrace = (h_shot[ix].imx)*nmy*nangx*nangy + (h_shot[ix].imy)*nangx*nangy
				position_m = 4*nz*itrace
				seek(stream_m,position_m)
				m = read(stream_m,Float32,nz*nangx*nangy)
				m = reshape(m,convert(Int64,nz),convert(Int64,nangx*nangy))
				for iz = 1 : nz
					# bilinear interpolation of angx and angy onto grid points
					angx = angx_shot[iz,ix]
					iangx = int(floor((angx - oangx)/dangx)) + 1
					angx1 = (iangx-1)*dangx + oangx
					angx2 = angx1 + dangx
					angy = angy_shot[iz,ix]
					iangy = int(floor((angy - oangy)/dangy)) + 1
					angy1 = (iangy-1)*dangy + oangy
					angy2 = angy1 + dangy
					w1 = (angx2-angx)*(angy2-angy)/((angx2-angx1)*(angy2-angy1))
					w2 = (angx-angx1)*(angy2-angy)/((angx2-angx1)*(angy2-angy1))
					w3 = (angx2-angx)*(angy-angy1)/((angx2-angx1)*(angy2-angy1))
					w4 = (angx-angx1)*(angy-angy1)/((angx2-angx1)*(angy2-angy1))
					if (iangx >= 1 && iangx <= nangx && iangy >= 1 && iangy <= nangy) 
						m[iz,(iangx-1)*nangy + iangy]   += w1*m_shot[iz,ix]
						if (iangx < nangx)
							m[iz,(iangx)*nangy   + iangy]   += w2*m_shot[iz,ix]
						end
						if (nangy > 1)
							m[iz,(iangx-1)*nangy + iangy+1] += w3*m_shot[iz,ix]
							if (iang < nang)
								m[iz,(iangx)*nangy   + iangy+1] += w4*m_shot[iz,ix]
							end
						end
					end
				end   
				seek(stream_m,position_m)
				write(stream_m,convert(Array{Float32,1},m[:]))
			end
			rm(join([shot_list[ishot].m ".seisd"]))
			rm(join([shot_list[ishot].m ".seish"]))
			rm(join([shot_list[ishot].d ".seisd"]))
			rm(join([shot_list[ishot].d ".seish"]))
			if (nangx != 1 || nangy != 1)
				rm(join([shot_list[ishot].angx ".seisd"]))
				rm(join([shot_list[ishot].angx ".seish"]))
				rm(join([shot_list[ishot].angy ".seisd"]))
				rm(join([shot_list[ishot].angy ".seish"]))
			end
		end	
		close(stream_m)
		close(stream_h)

	else
		m_m = join([m ".seisd"])
		m_h = join([m ".seish"])
		stream_m = open(m_m,"r")
		stream_h = open(m_h,"r")

		for ishot = 1 : nshot
			sx = shot_list[ishot].sx
			min_imxa = min_imx > int(floor((sx + ohx)/dhx)) - min_imx ? min_imx : int(floor((sx + ohx)/dhx)) - min_imx
			max_imxa = max_imx < int(floor((sx + ohx)/dhx)) + nhx - 1 ? max_imx : int(floor((sx + ohx)/dhx)) + nhx - 1
			nmxa = max_imxa - min_imxa + 1
			sy = shot_list[ishot].sy
			min_imya = min_imy > int(floor((sy + ohy)/dhy)) - min_imy ? min_imy : int(floor((sy + ohy)/dhy)) - min_imy
			max_imya = max_imy < int(floor((sy + ohy)/dhy)) + nhy - 1 ? max_imy : int(floor((sy + ohy)/dhy)) + nhy - 1
			nmya = max_imya - min_imya + 1
			j = 1    				
			h_shot = Array(Header,nmxa*nmya)
			for imx = 1 : nmxa
				for imy = 1 : nmya
					h_shot[j] = Seismic.InitSeisHeader() 
					h_shot[j].tracenum = convert(typeof(h_shot[1].tracenum),j)
					h_shot[j].o1 = convert(typeof(h_shot[1].o1),0)
					h_shot[j].n1 = convert(typeof(h_shot[1].n1),nz)
					h_shot[j].d1 = convert(typeof(h_shot[1].d1),dz)
					h_shot[j].imx = convert(typeof(h_shot[1].imx),imx-1 + min_imxa)
					h_shot[j].imy = convert(typeof(h_shot[1].imy),imy-1 + min_imya)
					h_shot[j].mx = convert(typeof(h_shot[1].mx),dmx*h_shot[j].imx)
					h_shot[j].my = convert(typeof(h_shot[1].my),dmy*h_shot[j].imy)
					j += 1
				end
			end
			m_shot = zeros(nz,nmxa*nmya)
			if (nangx != 1 || nangy != 1)
				angx_shot,h_ang = SeisRead(shot_list[ishot].angx)
				angy_shot,h_ang = SeisRead(shot_list[ishot].angy)
			else
				angx_shot = 0.*m_shot
				angy_shot  = 0.*m_shot
			end	
			for imx = 1 : nmxa
				for imy = 1 : nmya
					ix = (imx - 1)*nmy + imy
					h_shot[ix].sx = convert(typeof(h[1].sx),shot_list[ishot].sx)
					h_shot[ix].sy = convert(typeof(h[1].sy),shot_list[ishot].sy)	
					itrace = (imx-1 + min_imxa)*nmy*nangx*nangy + (imy-1 + min_imya)*nangx*nangy
					position_m = 4*nz*itrace
					seek(stream_m,position_m)
					m = read(stream_m,Float32,nz*nangx*nangy)
					m = reshape(m,convert(Int64,nz),convert(Int64,nangx*nangy))
					for iz = 1 : nz
						# bilinear interpolation of angx and angy onto grid points
						angx = angx_shot[iz,ix]
						iangx = int(floor((angx - oangx)/dangx)) + 1
						angx1 = (iangx-1)*dangx + oangx
						angx2 = angx1 + dangx
						angy = angy_shot[iz,ix]
						iangy = int(floor((angy - oangy)/dangy)) + 1
						angy1 = (iangy-1)*dangy + oangy
						angy2 = angy1 + dangy            		
						w1 = (angx2-angx)*(angy2-angy)/((angx2-angx1)*(angy2-angy1))
						w2 = (angx-angx1)*(angy2-angy)/((angx2-angx1)*(angy2-angy1))
						w3 = (angx2-angx)*(angy-angy1)/((angx2-angx1)*(angy2-angy1))
						w4 = (angx-angx1)*(angy-angy1)/((angx2-angx1)*(angy2-angy1))

						if (iangx >= 1 && iangx <= nangx && iangy >= 1 && iangy <= nangy) 
							m_shot[iz,ix] += w1*m[iz,(iangx-1)*nangy + iangy]
							if (iangx < nangx)
								m_shot[iz,ix] += w2*m[iz,(iangx)*nangy   + iangy]
							end
							if (nangy > 1)
								m_shot[iz,ix] += w3*m[iz,(iangx-1)*nangy + iangy+1]
								if (iangx < nangx)
									m_shot[iz,ix] += w4*m[iz,(iangx)*nangy   + iangy+1]
								end
							end
						end
					end   
				end
			end
			SeisWrite(shot_list[ishot].m,m_shot,h_shot)
		end	
		close(stream_m)
		close(stream_h)	
		a = pmap(shotwem,shot_list)
		j = 1
		for ishot = 1 : nshot
			d_shot,h_shot = SeisRead(shot_list[ishot].d)
			SeisWrite(d,d_shot,h_shot,["itrace"=>j])
			j += size(d_shot,2)
			SeisRemove(shot_list[ishot].d)
			SeisRemove(shot_list[ishot].m)
			if (nangx != 1 || nangy != 1)
				SeisRemove(shot_list[ishot].angx)
				SeisRemove(shot_list[ishot].angy)
			end			
		end	
	end
end

type Shot
	m
	d
	vel
	angx
	angy
	wav
	damping
	sx
	sy
	sz
	gz
	ohx
	dhx
	nhx
	ohy
	dhy
	nhy
	fmin
	fmax
	padt
	padx
	adj
	pade_flag
	verbose
end

function shotwem(shot)
	if (find_library(["shotwem"]) == "")
		error("Couldn't find shared library shotwem.so, make sure it is installed and check LD_LIBRARY_PATH environmental variable.")
	end
	argv = ["0",
	join(["adj=",shot.adj]), 
	join(["d=",shot.d]), join(["m=",shot.m]), join(["vel=",shot.vel]),  join(["wav=",shot.wav]),  
	join(["damping=",shot.damping]), join(["sx=",shot.sx]), join(["sy=",shot.sy]),  join(["sz=",shot.sz]),  join(["gz=",shot.gz]), 
	join(["ohx=",shot.ohx]),  join(["dhx=",shot.dhx]),  join(["nhx=",shot.nhx]), 
	join(["ohy=",shot.ohy]),  join(["dhy=",shot.dhy]),  join(["nhy=",shot.nhy]),   
	join(["fmin=",shot.fmin]),  join(["fmax=",shot.fmax]), 
	join(["padt=",shot.padt]),  join(["padx=",shot.padx]), 
	join(["pade_flag=",shot.pade_flag]),  join(["verbose=",shot.verbose]) ] 
	a = ccall((:main, "shotwem"), Int32, (Int32, Ptr{Ptr{Uint8}}), length(argv), argv)                    
	return(a)
end

