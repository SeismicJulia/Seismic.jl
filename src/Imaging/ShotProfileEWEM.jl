"""
**ShotProfileEWEM**

*Shot Profile Elastic Wave Equation Migration and Demigration of 3D isotropic data.*

**IN**   

* m : vector of filenames of image (comprised of [mpp;mps1;mps2])
* d : vector of filenames of data (comprised of [ux;uy;uz])
* adj : flag for adjoint (migration), or forward (demigration) (default=true) 
* damping = 1000. : damping for deconvolution imaging condition
* vp = "vp" : seis file containing the p-wave velocity (should have same x and z dimensions as the desired image)
* vs = "vs" : seis file containing the s-wave velocity (should have same x and z dimensions as the desired image)
* angx = "angx" : seis file containing incidence angles in the x direction for each shot
* angy = "angy" : seis file containing incidence angles in the y direction for each shot
* wav = "wav" : seis file containing the source wavelet (in time domain) 
* sz = 0. : source depth (Dev: read this from source wavelet file for variable source depth)
* gz = 0. : receiver depth (Dev: read this from data file for variable source depth (but then what to do in fwd op?))
* nhx = 101 : number of offset bins
* ohx = 1000. : min offset (surface offset in the data)
* dhx = 10. : offset increment
* nhy = 101 : number of offset bins
* ohy = 1000. : min offset (surface offset in the data)
* dhy = 10. : offset increment
* pade_flag = false : flag for Pade Fourier correction
* nangx = 1 : number of angle bins in x direction
* oangx = 0. : min angle in x direction (angle between source incidence angle and reflector normal in Degrees)
* dangx = 1. : angle increment in x direction
* nangy = 1 : number of angle bins in y direction
* oangy = 0. : min angle in y direction (angle between source incidence angle and reflector normal in Degrees)
* dangy = 1. : angle increment in y direction
* fmin = 0. : min frequency to process (Hz)
* fmax = 80. : max frequency to process (Hz)
* padt = 2 : pad factor for the time axis
* padx = 2 : pad factor for the spatial axes
* verbose = false : flag for error / debugging messages
* sx = [0.] : array of source X positions (meters)
* sy = [0.] : array of source Y positions (meters)

**OUT**  

*Credits: AS, 2015*

"""

function ShotProfileEWEM(m::Array{ASCIIString,1},d::Array{ASCIIString,1},adj=true;damping=1000.,vp="vp.seis",vs="vs.seis",angx="angx.seis",angy="angy.seis",wav="wav.seis",sz=0.,gz=0.,nhx=100,ohx=0,dhx=10,nhy=1,ohy=0,dhy=10,pade_flag=false,nangx=1,oangx=0,dangx=1,nangy=1,oangy=0,dangy=1,fmin=0,fmax=80,padt=2,padx=2,verbose=false,sx=[0],sy=[0])


	nshot = length(sx)	
	vel,h,extent = SeisRead(vp)
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

	shot_list = Array(Shot3C,nshot)
	for ishot = 1 : nshot
		shot_list[ishot] = Shot3C(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		shot_list[ishot].ux = join([d[1] "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].uy = join([d[2] "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].uz = join([d[3] "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].mpp = join([m[1] "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].mps1 = join([m[2] "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].mps2 = join([m[3] "_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].angx = join(["angx_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].angy = join(["angy_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].vp = vp
		shot_list[ishot].vs = vs
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
			SeisWindow(angx,shot_list[ishot].angx,key=["sx","sy"],minval=[sx[ishot],sy[ishot]],maxval=[sx[ishot],sy[ishot]])
			SeisWindow(angy,shot_list[ishot].angy,key=["sx","sy"],minval=[sx[ishot],sy[ishot]],maxval=[sx[ishot],sy[ishot]])
		end
	end

	if (adj == true)
		for ishot = 1 : nshot
			SeisWindow(d[1],shot_list[ishot].ux,key=["sx","sy","gx","gy"],minval=[sx[ishot],sy[ishot],min_gx,min_gy],maxval=[sx[ishot],sy[ishot],max_gx,max_gy])
			SeisWindow(d[2],shot_list[ishot].uy,key=["sx","sy","gx","gy"],minval=[sx[ishot],sy[ishot],min_gx,min_gy],maxval=[sx[ishot],sy[ishot],max_gx,max_gy])
			SeisWindow(d[3],shot_list[ishot].uz,key=["sx","sy","gx","gy"],minval=[sx[ishot],sy[ishot],min_gx,min_gy],maxval=[sx[ishot],sy[ishot],max_gx,max_gy])
		end
		a = pmap(shotewem,shot_list)
		j = 1
		gather = zeros(Float32,nz,nangx*nangy)
		extent = Seismic.Extent(nz,nmx,nmy,nangx,nangy,
				oz,dhx*min_imx,dhy*min_imy,oangx,oangy,
				dz,dhx,dhy,dangx,dangy,
				"Depth","mx","my","angx","angy",
				"s","m","m","degrees","degrees",
				"")
		for imx = 1 : nmx
			for imy = 1 : nmy
				h = Header[]
				for iangx = 1 : nangx
					for iangy = 1 : nangy
						push!(h,Seismic.InitSeisHeader())		
						h[(iangx-1)*nangy + iangy].tracenum = convert(typeof(h[1].tracenum),j + iangx*(nangy-1) + iangy)
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
				SeisWrite(m[1],gather,h,extent,itrace=j)
				SeisWrite(m[2],gather,h,extent,itrace=j)
				SeisWrite(m[3],gather,h,extent,itrace=j)
				j += nangx*nangy
			end
		end

		mpp_m = ParseDataName(m[1])
		mpp_h = ParseHeaderName(m[1])
		mps1_m = ParseDataName(m[2])
		mps1_h = ParseHeaderName(m[2])
		mps2_m = ParseDataName(m[3])
		mps2_h = ParseHeaderName(m[3])
		stream_mpp = open(mpp_m,"a+")
		stream_hpp = open(mpp_h,"a+")
		stream_mps1 = open(mps1_m,"a+")
		stream_hps1 = open(mps1_h,"a+")
		stream_mps2 = open(mps2_m,"a+")
		stream_hps2 = open(mps2_h,"a+")

		for ishot = 1 : nshot
			mpp_shot,hpp_shot,extent  = SeisRead(shot_list[ishot].mpp)
			mps1_shot,hps1_shot,extent  = SeisRead(shot_list[ishot].mps1)
			mps2_shot,hps2_shot,extent  = SeisRead(shot_list[ishot].mps2)
			if (nangx != 1 || nangy != 1)
				angx_shot,h_ang,extent = SeisRead(shot_list[ishot].angx)
				angy_shot,h_ang,extent = SeisRead(shot_list[ishot].angy)
			else
				angx_shot = 0.*mpp_shot
				angy_shot = 0.*mpp_shot
			end	

			nx_shot = size(mpp_shot,2)
			for ix = 1 : nx_shot
				itrace = (hpp_shot[ix].imx)*nmy*nangx*nangy + (hpp_shot[ix].imy)*nangx*nangy
				position_m = 4*nz*itrace
				seek(stream_mpp,position_m)
				m1 = read(stream_mpp,Float32,nz*nangx*nangy)
				m1 = reshape(m1,convert(Int64,nz),convert(Int64,nangx*nangy)) 
				seek(stream_mps1,position_m)
				m2 = read(stream_mps1,Float32,nz*nangx*nangy)
				m2 = reshape(m2,convert(Int64,nz),convert(Int64,nangx*nangy)) 
				seek(stream_mps2,position_m)
				m3 = read(stream_mps2,Float32,nz*nangx*nangy)
				m3 = reshape(m3,convert(Int64,nz),convert(Int64,nangx*nangy)) 
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
						m1[iz,(iangx-1)*nangy + iangy]   += w1*mpp_shot[iz,ix]
						if (iangx < nangx)
							m1[iz,(iangx)*nangy   + iangy]   += w2*mpp_shot[iz,ix]
						end
						if (nangy > 1)
							m1[iz,(iangx-1)*nangy + iangy+1] += w3*mpp_shot[iz,ix]
							if (iangx < nangx)
								m1[iz,(iangx)*nangy   + iangy+1] += w4*mpp_shot[iz,ix]
							end
						end
						m2[iz,(iangx-1)*nangy + iangy]   += w1*signf(angx)*signf(angy)*mps1_shot[iz,ix]
						if (iangx < nangx)
							m2[iz,(iangx)*nangy   + iangy]   += w2*signf(angx)*signf(angy)*mps1_shot[iz,ix]
						end
						if (nangy > 1)
							m2[iz,(iangx-1)*nangy + iangy+1] += w3*signf(angx)*signf(angy)*mps1_shot[iz,ix]
							if (iangx < nangx)
								m2[iz,(iangx)*nangy   + iangy+1] += w4*signf(angx)*signf(angy)*mps1_shot[iz,ix]
							end
						end
						m3[iz,(iangx-1)*nangy + iangy]   += w1*signf(angx)*signf(angy)*mps2_shot[iz,ix]
						if (iangx < nangx)
							m3[iz,(iangx)*nangy   + iangy]   += w2*signf(angx)*signf(angy)*mps2_shot[iz,ix]
						end
						if (nangy > 1)
							m3[iz,(iangx-1)*nangy + iangy+1] += w3*signf(angx)*signf(angy)*mps2_shot[iz,ix]
							if (iangx < nangx)
								m3[iz,(iangx)*nangy   + iangy+1] += w4*signf(angx)*signf(angy)*mps2_shot[iz,ix]
							end
						end
					end
				end   
				seek(stream_mpp,position_m)
				write(stream_mpp,convert(Array{Float32,1},m1[:]))
				seek(stream_mps1,position_m)
				write(stream_mps1,convert(Array{Float32,1},m2[:]))
				seek(stream_mps2,position_m)
				write(stream_mps2,convert(Array{Float32,1},m3[:]))
			end
			SeisRemove(shot_list[ishot].ux)
			SeisRemove(shot_list[ishot].uy)
			SeisRemove(shot_list[ishot].uz)
			SeisRemove(shot_list[ishot].mpp)
			SeisRemove(shot_list[ishot].mps1)
			SeisRemove(shot_list[ishot].mps2)			
			if (nangx != 1 || nangy != 1)
				SeisRemove(shot_list[ishot].angx)
				SeisRemove(shot_list[ishot].angy)
			end
		end	
		close(stream_mpp)
		close(stream_hpp)
		close(stream_mps1)
		close(stream_hps1)
		close(stream_mps2)
		close(stream_hps2)
	else

		mpp_m = join([m[1] ".seisd"])
		mpp_h = join([m[1] ".seish"])
		mps1_m = join([m[2] ".seisd"])
		mps1_h = join([m[2] ".seish"])
		mps2_m = join([m[3] ".seisd"])
		mps2_h = join([m[3] ".seish"])
		stream_mpp = open(mpp_m,"r")
		stream_hpp = open(mpp_h,"r")
		stream_mps1 = open(mps1_m,"r")
		stream_hps1 = open(mps1_h,"r")
		stream_mps2 = open(mps2_m,"r")
		stream_hps2 = open(mps2_h,"r")

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

			mpp_shot = zeros(nz,nmxa*nmya)
			mps1_shot = zeros(nz,nmxa*nmya)
			mps2_shot = zeros(nz,nmxa*nmya)
			if (nangx != 1 || nangy != 1)
				angx_shot,h_ang = SeisRead(shot_list[ishot].angx)
				angy_shot,h_ang = SeisRead(shot_list[ishot].angy)
			else
				angx_shot = 0.*mpp_shot
				angy_shot = 0.*mpp_shot
			end	
			for imx = 1 : nmxa
				for imy = 1 : nmya
					ix = (imx - 1)*nmya + imy
					h_shot[ix].sx = convert(typeof(h[1].sx),shot_list[ishot].sx)
					h_shot[ix].sy = convert(typeof(h[1].sy),shot_list[ishot].sy)	
					itrace = (imx-1 + min_imxa)*nmy*nangx*nangy + (imy-1 + min_imya)*nangx*nangy
					position_m = 4*nz*itrace
					seek(stream_mpp,position_m)
					m1 = read(stream_mpp,Float32,nz*nangx*nangy)
					m1 = reshape(m1,convert(Int64,nz),convert(Int64,nangx*nangy))
					seek(stream_mps1,position_m)
					m2 = read(stream_mps1,Float32,nz*nangx*nangy)
					m2 = reshape(m2,convert(Int64,nz),convert(Int64,nangx*nangy))
					seek(stream_mps2,position_m)
					m3 = read(stream_mps2,Float32,nz*nangx*nangy)
					m3 = reshape(m3,convert(Int64,nz),convert(Int64,nangx*nangy))
					for iz = 1 : nz
						# bilinear interpolation of ang and az onto grid points
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
							mpp_shot[iz,ix] += w1*m1[iz,(iangx-1)*nangy + iangy]
							if (iangx < nangx)
								mpp_shot[iz,ix] += w2*m1[iz,(iangx)*nangy   + iangy]
							end
							if (nangy > 1)
								mpp_shot[iz,ix] += w3*m1[iz,(iangx-1)*nangy + iangy+1]
								if (iangx < nangx)
									mpp_shot[iz,ix] += w4*m1[iz,(iangx)*nangy   + iangy+1]
								end
							end
							mps1_shot[iz,ix] += w1*signf(angx)*signf(angy)*m2[iz,(iangx-1)*nangy + iangy]
							if (iangx < nangx)
								mps1_shot[iz,ix] += w2*signf(angx)*signf(angy)*m2[iz,(iangx)*nangy   + iangy]
							end
							if (nangy > 1)
								mps1_shot[iz,ix] += w3*signf(angx)*signf(angy)*m2[iz,(iangx-1)*nangy + iangy+1]
								if (iangx < nangx)
									mps1_shot[iz,ix] += w4*signf(angx)*signf(angy)*m2[iz,(iangx)*nangy   + iangy+1]
								end
							end
							mps2_shot[iz,ix] += w1*signf(angx)*signf(angy)*m3[iz,(iangx-1)*nangy + iangy]
							if (iangx < nangx)
								mps2_shot[iz,ix] += w2*signf(angx)*signf(angy)*m3[iz,(iangx)*nangy   + iangy]
							end
							if (nangy > 1)
								mps2_shot[iz,ix] += w3*signf(angx)*signf(angy)*m3[iz,(iangx-1)*nangy + iangy+1]
								if (iangx < nangx)
									mps2_shot[iz,ix] += w4*signf(angx)*signf(angy)*m3[iz,(iangx)*nangy   + iangy+1]
								end
							end
						end
					end   
				end
			end
			SeisWrite(shot_list[ishot].mpp,mpp_shot,h_shot)
			SeisWrite(shot_list[ishot].mps1,mps1_shot,h_shot)
			SeisWrite(shot_list[ishot].mps2,mps2_shot,h_shot)
		end	
		close(stream_mpp)
		close(stream_hpp)	
		close(stream_mps1)
		close(stream_hps1)	
		close(stream_mps2)
		close(stream_hps2)	
		a = pmap(shotewem,shot_list)
		j = 1
		for ishot = 1 : nshot
			ux_shot,h_shot = SeisRead(shot_list[ishot].ux)
			SeisWrite(d[1],ux_shot,h_shot,itrace=j)
			uy_shot,h_shot = SeisRead(shot_list[ishot].uy)
			SeisWrite(d[2],uy_shot,h_shot,itrace=j)
			uz_shot,h_shot = SeisRead(shot_list[ishot].uz)
			SeisWrite(d[3],uz_shot,h_shot,itrace=j)
			j += size(ux_shot,2)
			SeisRemove(shot_list[ishot].ux)
			SeisRemove(shot_list[ishot].uy)
			SeisRemove(shot_list[ishot].uz)
			SeisRemove(shot_list[ishot].mpp)
			SeisRemove(shot_list[ishot].mps1)
			SeisRemove(shot_list[ishot].mps2)
			if (nangx != 1 || nangy != 1)
				SeisRemove(shot_list[ishot].angx)
				SeisRemove(shot_list[ishot].angy)
			end
		end	
	end
end

type Shot3C
	ux
	uy
	uz
	mpp
	mps1
	mps2
	angx
	angy
	vp
	vs
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

function shotewem(shot)
	if (find_library(["shotewem"]) == "")
		error("Couldn't find shared library shotewem.so, make sure it is installed and check LD_LIBRARY_PATH environmental variable.")
	end
	argv = ["0",
	join(["adj=",shot.adj]), 
	join(["ux=",shot.ux]), join(["uy=",shot.uy]), join(["uz=",shot.uz]),
	join(["mpp=",shot.mpp]), join(["mps1=",shot.mps1]), join(["mps2=",shot.mps2]),
	join(["vp=",shot.vp]), join(["vs=",shot.vs]), join(["wav=",shot.wav]), 
	join(["damping=",shot.damping]), join(["sx=",shot.sx]), join(["sy=",shot.sy]),  join(["sz=",shot.sz]),  join(["gz=",shot.gz]), 
	join(["ohx=",shot.ohx]),  join(["dhx=",shot.dhx]),  join(["nhx=",shot.nhx]), 
	join(["ohy=",shot.ohy]),  join(["dhy=",shot.dhy]),  join(["nhy=",shot.nhy]),   
	join(["fmin=",shot.fmin]),  join(["fmax=",shot.fmax]), 
	join(["padt=",shot.padt]),  join(["padx=",shot.padx]), 
	join(["pade_flag=",shot.pade_flag]),  join(["verbose=",shot.verbose]) ] 
	a = ccall((:main, "shotewem"), Int32, (Int32, Ptr{Ptr{Uint8}}), length(argv), argv)                    
	return(a)

end

function signf(a)

	b = a < 0.0 ? -1.0 : 1.0
	return(b)

end

