function smooth_angles(d,h,param)
	
	Nsmooth = get(param,"Nsmooth",3)
	Nrepeat = get(param,"Nrepeat",1)
	min_iang = h[1].iang
	max_iang = h[end].iang
	min_iaz = h[1].iaz
	max_iaz = h[end].iaz
	nang = max_iang - min_iang + 1
	naz = max_iaz - min_iaz + 1
	nz = convert(Int64,h[1].n1)
	a = vec(zeros(nang,1))
	d = reshape(d,nz,nang,naz)	
	for iz = 1 : nz
		for iaz = 1 : naz
			for iang = 1 : nang
				a[iang] = d[iz,iang,iaz]
			end
			for irepeat = 1 : Nrepeat
				a = mean_filter(a,nang,Nsmooth)
			end
			for iang = 1 : nang
				d[iz,iang,iaz] = a[iang]
			end
		end
	end
	
	return d[1:nz,:],h;
end

function mean_filter(a,nx,nxw)

	b = vec(zeros(nx,1))
	for ix = 1 : nx
		sum  = 0.
		nsum = 0.
		for ixw = 1 : nxw
			index1 = ix - floor(nxw/2) - 1 + ixw
			if (index1>0 && index1<nx)
				sum += a[index1]
				nsum += 1.
			end
		end
    	b[ix] = sum/nsum
    end
	return b
end

function SmoothGathers(m::ASCIIString,d::ASCIIString,param=Dict())
	if (param["adj"]==false)
		SeisProcess(m,d,[smooth_angles],param,"gather",["imx","imy"])
	else
		SeisProcess(d,m,[smooth_angles],param,"gather",["imx","imy"])
	end
end

