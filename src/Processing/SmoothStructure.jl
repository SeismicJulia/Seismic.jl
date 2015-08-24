function smooth_structure(d,h,param)

	Nsmooth2 = get(param,"Nsmooth2",3)
	dz = 1
	dx = 1
	dip_field = get(param,"dip","dip")
	dip,h_dip,s = SeisRead(dip_field)
	dip2 = atan(abs(dip))

	nx = size(d,2)
	nz = convert(Int64,h[1].n1)
	WL = Nsmooth2*dx/2
	d2 = copy(d)	
	for iz = 1 : nz
		for ix = 1 : nx
			a = vec(zeros(Nsmooth2,1))
			#println("dip2=",dip2[iz,ix]*180/pi)
			ixo = ix - convert(Int32,floor(WL*cos(dip2[iz,ix])/dx))
			dip[iz,ix] > 0 ? izo = iz - convert(Int32,floor(WL*sin(dip2[iz,ix])/dz)) : izo = iz + convert(Int32,floor(WL*sin(dip2[iz,ix])/dz))
			for iw = 1 : Nsmooth2
				ixw = ixo + convert(Int32,floor(iw*dx*cos(dip2[iz,ix])/dx))
				dip[iz,ix] > 0 ? izw = izo + convert(Int32,floor(iw*dx*sin(dip2[iz,ix])/dz)) : izw = izo - convert(Int32,floor(iw*dx*sin(dip2[iz,ix])/dz))
				if (izw > 0 && izw < nz && ixw > 0 && ixw < nx)
					a[iw] = d[izw,ixw]
				end
			end
			d2[iz,ix] = sum(a[:])/Nsmooth2
		end
	end

	return d2[1:nz,:],h;
end

function SmoothStructure(m::ASCIIString,d::ASCIIString,param=Dict())
        if (param["adj"]==false)
		SeisSort(m,"tmp_smoothstruct1",["iang","imx"],false);
		SeisProcess("tmp_smoothstruct1","tmp_smoothstruct2",[smooth_structure],param,"gather",["iang"])
		SeisSort("tmp_smoothstruct2",d,["imx","iang"],false);	
	else
		SeisSort(d,"tmp_smoothstruct1",["iang","imx"],false);
		SeisProcess("tmp_smoothstruct1","tmp_smoothstruct2",[smooth_structure],param,"gather",["iang"])
		SeisSort("tmp_smoothstruct2",m,["imx","iang"],false);	
	end
	rm(join(["tmp_smoothstruct1" ".seisd"]))
	rm(join(["tmp_smoothstruct1" ".seish"]))
	rm(join(["tmp_smoothstruct2" ".seisd"]))
	rm(join(["tmp_smoothstruct2" ".seish"]))
end
