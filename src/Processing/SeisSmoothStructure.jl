function SeisSmoothStructure(m::AbstractString,d::AbstractString,param)

	adj = get(param,"adj",false)
	param["f"] = [smooth_structure]
	param["group"] = "gather"
	param["key"] = ["iaz","iang"]
	tmp1 = join(["tmp_SmoothStructure1_",string(int(rand()*100000))])
	tmp2 = join(["tmp_SmoothStructure2_",string(int(rand()*100000))])
	if (adj==false)
		SeisSort(m,tmp1,key=["iaz","iang","imy","imx"]);
		SeisProcess(tmp1,tmp2,param)
		SeisSort(tmp2,d,key=["imx","imy","iang","iaz"]);
	else
		SeisSort(d,tmp1,key=["iaz","iang","imy","imx"]);
		SeisProcess(tmp1,tmp2,param)
		SeisSort(tmp2,m,key=["imx","imy","iang","iaz"]);
	end
	SeisRemove(tmp1)
	SeisRemove(tmp2)

end

function smooth_structure(d,h,param)

	Nsmooth2 = get(param,"Nsmooth2",3)
	dip_field = get(param,"dip","dip")
	dip,h_dip = SeisRead(dip_field)
	#dip2 = atan(abs(dip))
	dip2 = atan(dip)

	nx = size(d,2)
	nz = size(d,1)
	WL = int(Nsmooth2/2)
	d2 = copy(d)
	smootharray!(d,d2,dip,dip2,nx,nz,Nsmooth2,WL)

	return d2,h;
end

function smootharray!(d::Array{Float32,2},d2::Array{Float32,2},dip::Array{Float32,2},dip2::Array{Float32,2},nx::Int64,nz::Int64,Nsmooth2::Int64,WL::Int64)

	for ix = WL : size(d,2) - WL
		for iz = WL : size(d,1) - WL
			@inbounds d2[iz,ix] = 0.
			@inbounds d2[iz,ix] = smooth1pixel(d,dip[iz,ix],cos(dip2[iz,ix]),sin(dip2[iz,ix]),ix - int(WL*cos(abs(dip2[iz,ix]))),iz - int(WL*sin(dip2[iz,ix])),nx,nz,Nsmooth2)
		end
	end

end

function smooth1pixel(m::Array{Float32,2},dip::Float32,xval::Float32,zval::Float32,ixo::Int64,izo::Int64,nx::Int64,nz::Int64,nw::Int64)

	val = 0.
	for iw = 1 : nw
		@inbounds val += m[izo + int(iw*zval),ixo + int(iw*xval)]
	end

	return val/nw

end
