"""
**ComputeAngles**

*Compute angles for shot gathers. These angles can be used for mapping migrated shots into angle gathers during shot profile migration.*

**IN**


* angx = "angx" : filename for incidence angles in the x direction for each shot
* angy = "angy" : filename for incidence angles in the y direction for each shot
* dip_flag = false : flag to subtract reflector dip from the computed angles to make them with reference to reflector normal
* vel = "vel" : seis file containing the velocity (should have same x and z dimensions as the desired image)
* wav = "wav" : seis file containing the source wavelet (in time domain)
* sz = 0. : source depth (Dev: read this from source wavelet file for variable source depth)
* nhx = 101 : number of offset bins
* ohx = 1000. : min offset (surface offset in the data)
* dhx = 10. : offset increment
* nhy = 101 : number of offset bins
* ohy = 1000. : min offset (surface offset in the data)
* dhy = 10. : offset increment
* pade_flag = false : flag for Pade Fourier correction
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

function ComputeAngles(angx::AbstractString,angy::AbstractString,dip_flag=false,vel="vel",wav="wav",sz=0.,nhx=100,ohx=0,dhx=10,nhy=1,ohy=0,dhy=10,pade_flag=false,fmin=0,fmax=80,padt=2,padx=2,verbose=false,sx=[0],sy=[0])

	nshot = length(sx)
	if (dip_flag=="y")
		dipx,h = SeisRead(dipx_name)
		dipy,h = SeisRead(dipy_name)
	end

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

	shot_list = Array(ShotAngles,nshot)
	for ishot = 1 : nshot
		shot_list[ishot] = ShotAngles(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		shot_list[ishot].angx = join(["angx_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].angy = join(["angy_shot_" int(floor(sx[ishot])) "_" int(floor(sy[ishot]))])
		shot_list[ishot].vel = vel
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
		shot_list[ishot].dip_flag = (dip_flag == true) ? "y" : "n"
		shot_list[ishot].verbose = (verbose == true) ? "y" : "n"
	end

	a = pmap(compute_angles,shot_list)

	j = 1
	for ishot = 1 : nshot
		angx_shot,h_shot = SeisRead(shot_list[ishot].angx)
		angy_shot,h_shot = SeisRead(shot_list[ishot].angy)
		SeisWrite(angx,angx_shot,h_shot,itrace=j)
		SeisWrite(angy,angy_shot,h_shot,itrace=j)
		j += size(angx_shot,2)
		SeisRemove(shot_list[ishot].angx)
		SeisRemove(shot_list[ishot].angy)
	end

end

type ShotAngles
	angx
	angy
	vel
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

function compute_angles(shot)
	ENV["OMP_NUM_THREADS"] = shot.omp
	if (find_library(["compute_angles"]) == "")
		error("Couldn't find shared library compute_angles.so, make sure it is installed and check LD_LIBRARY_PATH environmental variable.")
	end
	argv = ["0",
	join(["angx=",shot.angx]), join(["angy=",shot.angy]), join(["vel=",shot.vel]), join(["dipx=",shot.dipx]), join(["dipy=",shot.dipy]), join(["wav=",shot.wav]),
	join(["sx=",shot.sx]), join(["sy=",shot.sy]),  join(["sz=",shot.sz]),
	join(["ohx=",shot.ohx]),  join(["dhx=",shot.dhx]),  join(["nhx=",shot.nhx]),
	join(["ohy=",shot.ohy]),  join(["dhy=",shot.dhy]),  join(["nhy=",shot.nhy]),
	join(["fmin=",shot.fmin]),  join(["fmax=",shot.fmax]),
	join(["pade_flag=",shot.pade_flag]),  join(["dip_flag=",shot.dip_flag]),  join(["verbose=",shot.verbose]) ]
	a = ccall((:main, "compute_angles"), Int32, (Int32, Ptr{Ptr{UInt8}}), length(argv), argv)
	return(a)

end
