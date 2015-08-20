include("Header.jl")

function Taper(d,nt,nx1,nx2,nx3,nx4,tapti,taptf,tapx1i,tapx1f,tapx2i,tapx2f,tapx3i,tapx3f,tapx4i,tapx4f)


	#println("nt=",nt)
	#println("tapti=",tapti)
	#println("taptf=",taptf)

	tx1=1;tx2=1;tx3=1;tx4=1;
	for ix1 = 1 : nx1
    	if (ix1>=1   && ix1<tapx1i) 
    		tx1 = 1 - cos(((ix1-1)/tapx1i)*pi/2)
	end
    	if (ix1>=tapx1i && ix1<=nx1-tapx1f) 
    		tx1 = 1
	end
    	if (ix1>nx1-tapx1f && ix1<=nx1)
    		tx1 = cos(((ix1-1-nx1+tapx1f)/tapx1f)*pi/2)
    	end
  		for ix2 = 1 : nx2
    		if (ix2>=1   && ix2<tapx2i) 
    			tx2 = 1 - cos(((ix2-1)/tapx2i)*pi/2)
		end
    		if (ix2>=tapx2i && ix2<=nx2-tapx2f) 
    			tx2 = 1
		end
    		if (ix2>nx2-tapx2f && ix2<=nx2)
    			tx2 = cos(((ix2-1-nx2+tapx2f)/tapx2f)*pi/2)
    		end
  			for ix3 = 1 : nx3
				if (ix3>=1   && ix3<tapx3i) 
    					tx3 = 1 - cos(((ix3-1)/tapx3i)*pi/2)
				end
    			        if (ix3>=tapx3i && ix3<=nx3-tapx3f) 
    					tx3 = 1
				end
    			        if (ix3>nx3-tapx3f && ix3<=nx3)
    					tx3 = cos(((ix3-1-nx3+tapx3f)/tapx3f)*pi/2)
    				end
  				for ix4 = 1 : nx4
    				if (ix4>=1   && ix4<tapx4i) 
    					tx4 = 1 - cos(((ix4-1)/tapx4i)*pi/2)
				end
    				if (ix4>=tapx4i && ix4<=nx4-tapx4f)
    					tx4 = 1
				end
    				if (ix4>nx4-tapx4f && ix4<=nx4)
    					tx4 = cos(((ix4-1-nx4+tapx4f)/tapx4f)*pi/2)
    				end	
    				ix = (ix1-1)*nx2*nx3*nx4 + (ix2-1)*nx3*nx4 + (ix3-1)*nx4 + ix4
  					for it = 1 : nt
    					if (it>=1 && it<tapti) 
    						tt = 1 - cos(((it-1)/tapti)*pi/2)
					end		
    					if (it>=tapti && it<=nt-taptf) 
    						tt = 1
					end
    					if (it>nt-taptf && it<=nt)
    						tt = cos(((it-1-nt+taptf)/taptf)*pi/2)
    					end
    					d[it,ix] = tt*tx1*tx2*tx3*tx4*d[it,ix]
    					#if (ix==1) 
    					#	println("tt ",it,"= ",tt)
    					#end
  					end
  				end
  			end
  		end
	end
	return d
end

function SeisUnPatch(list,out,param)

    style = get(param,"style","sxsygxgy")
    if (style == "sxsygxgy")
		key = ["t","isx","isy","igx","igy"]
		min_ix1 = get(param,"min_isx",0)
		max_ix1 = get(param,"max_isx",0)
		min_ix2 = get(param,"min_isy",0)
		max_ix2 = get(param,"max_isy",0)
		min_ix3 = get(param,"min_igx",0)
		max_ix3 = get(param,"max_igx",0)
		min_ix4 = get(param,"min_igy",0)
		max_ix4 = get(param,"max_igy",0)
	elseif (style=="mxmyhxhy")
		key = ["t","imx","imy","ihx","ihy"]
		min_ix1 = get(param,"min_imx",0)
		max_ix1 = get(param,"max_imx",0)
		min_ix2 = get(param,"min_imy",0)
		max_ix2 = get(param,"max_imy",0)
		min_ix3 = get(param,"min_ihx",0)
		max_ix3 = get(param,"max_ihx",0)
		min_ix4 = get(param,"min_ihy",0)
		max_ix4 = get(param,"max_ihy",0)
	elseif (style=="mxmyhaz")
		key = ["t","imx","imy","ih","iaz"]
		min_ix1 = get(param,"min_imx",0)
		max_ix1 = get(param,"max_imx",0)
		min_ix2 = get(param,"min_imy",0)
		max_ix2 = get(param,"max_imy",0)
		min_ix3 = get(param,"min_ih",0)
		max_ix3 = get(param,"max_ih",0)
		min_ix4 = get(param,"min_iaz",0)
		max_ix4 = get(param,"max_iaz",0)
	elseif (style=="sxsyhxhy")
		key = ["t","isx","isy","ihx","ihy"]
		min_ix1 = get(param,"min_isx",0)
		max_ix1 = get(param,"max_isx",0)
		min_ix2 = get(param,"min_isy",0)
		max_ix2 = get(param,"max_isy",0)
		min_ix3 = get(param,"min_ihx",0)
		max_ix3 = get(param,"max_ihx",0)
		min_ix4 = get(param,"min_ihy",0)
		max_ix4 = get(param,"max_ihy",0)
	elseif (style=="gxgyhxhy")
		key = ["t","igx","igy","ihx","ihy"]
		min_ix1 = get(param,"min_igx",0)
		max_ix1 = get(param,"max_igx",0)
		min_ix2 = get(param,"min_igy",0)
		max_ix2 = get(param,"max_igy",0)
		min_ix3 = get(param,"min_ihx",0)
		max_ix3 = get(param,"max_ihx",0)
		min_ix4 = get(param,"min_ihy",0)
		max_ix4 = get(param,"max_ihy",0)
	elseif (style=="sxsyhaz")
		key = ["t","isx","isy","ih","iaz"]
		min_ix1 = get(param,"min_isx",0)
		max_ix1 = get(param,"max_isx",0)
		min_ix2 = get(param,"min_isy",0)
		max_ix2 = get(param,"max_isy",0)
		min_ix3 = get(param,"min_ih",0)
		max_ix3 = get(param,"max_ih",0)
		min_ix4 = get(param,"min_iaz",0)
		max_ix4 = get(param,"max_iaz",0)
	elseif (style=="gxgyhaz")
		key = ["t","igx","igy","ih","iaz"]
		min_ix1 = get(param,"min_igx",0)
		max_ix1 = get(param,"max_igx",0)
		min_ix2 = get(param,"min_igy",0)
		max_ix2 = get(param,"max_igy",0)
		min_ix3 = get(param,"min_ih",0)
		max_ix3 = get(param,"max_ih",0)
		min_ix4 = get(param,"min_iaz",0)
		max_ix4 = get(param,"max_iaz",0)
	else
		error("style not defined.")
	end
	nx1 = max_ix1 - min_ix1 + 1
	nx2 = max_ix2 - min_ix2 + 1
	nx3 = max_ix3 - min_ix3 + 1
	nx4 = max_ix4 - min_ix4 + 1
	
	nt = get(param,"nt",0)
	it_WL = get(param,"it_WL",nt)
	it_WO = get(param,"it_WO",0)
	ix1_WL = get(param,"ix1_WL",nx1)
	ix1_WO = get(param,"ix1_WO",0)
	ix2_WL = get(param,"ix2_WL",nx2)
	ix2_WO = get(param,"ix2_WO",0)
	ix3_WL = get(param,"ix3_WL",nx3)
	ix3_WO = get(param,"ix3_WO",0)
	ix4_WL = get(param,"ix4_WL",nx4)
	ix4_WO = get(param,"ix4_WO",0)

	
	ang = get(param,"ang",90)
	gamma = get(param,"gamma",1)
	rad2deg = 180/pi;
	deg2rad = pi/180;
        gammainv = 1/gamma;
        if (ang > 90) 
    	    ang2=-deg2rad*(ang-90)
        else 
    	    ang2=deg2rad*(90-ang)
        end

	osx = get(param,"osx",0)
	osy = get(param,"osy",0)
	ogx = get(param,"ogx",0)
	ogy = get(param,"ogy",0)
	omx = get(param,"omx",0)
	omy = get(param,"omy",0)
	ohx = get(param,"ohx",0)
	ohy = get(param,"ohy",0)
	oh = get(param,"oh",0)
	oaz = get(param,"oaz",0)

	dsx = get(param,"dsx",1)
	dsy = get(param,"dsy",1)
	dgx = get(param,"dgx",1)
	dgy = get(param,"dgy",1)
	dmx = get(param,"dmx",1)
	dmy = get(param,"dmy",1)
	dhx = get(param,"dhx",1)
	dhy = get(param,"dhy",1)
	dh  = get(param,"dh",1)
	daz = get(param,"daz",1)

	min_isx = get(param,"min_isx",0)
	max_isx = get(param,"max_isx",0)
    nsx = max_isx - min_isx + 1
	min_isy = get(param,"min_isy",0)
	max_isy = get(param,"max_isy",0)
    nsy = max_isy - min_isy + 1
	min_igx = get(param,"min_igx",0)
	max_igx = get(param,"max_igx",0)
    ngx = max_igx - min_igx + 1
	min_igy = get(param,"min_igy",0)
	max_igy = get(param,"max_igy",0)
    ngy = max_igy - min_igy + 1

	min_imx = get(param,"min_imx",0)
	max_imx = get(param,"max_imx",0)
    nmx = max_imx - min_imx + 1
	min_imy = get(param,"min_imy",0)
	max_imy = get(param,"max_imy",0)
    nmy = max_imy - min_imy + 1
	min_ihx = get(param,"min_ihx",0)
	max_ihx = get(param,"max_ihx",0)
    nhx = max_ihx - min_ihx + 1
	min_ihy = get(param,"min_ihy",0)
	max_ihy = get(param,"max_ihy",0)
    nhy = max_ihy - min_ihy + 1
	min_ih = get(param,"min_ih",0)
	max_ih = get(param,"max_ih",0)
    nh = max_ih - min_ih + 1
	min_iaz = get(param,"min_iaz",0)
	max_iaz = get(param,"max_iaz",0)
    naz = max_iaz - min_iaz + 1

    if (style=="sxsygxgy")
	nx1=nsx;nx2=nsy;nx3=ngx;nx4=ngy
    elseif (style=="mxmyhxhy")
	nx1=nmx;nx2=nmy;nx3=nhx;nx4=nhy
    elseif (style=="mxmyhaz")
	nx1=nmx;nx2=nmy;nx3=nh;nx4=naz
    elseif (style=="sxsyhxhy")
	nx1=nsx;nx2=nsy;nx3=nhx;nx4=nhy
    elseif (style=="gxgyhxhy")
	nx1=ngx;nx2=ngy;nx3=nhx;nx4=nhy
    elseif (style=="sxsyhaz")
	nx1=nsx;nx2=nsy;nx3=nh;nx4=naz
    elseif (style=="gxgyhaz")
	nx1=ngx;nx2=ngy;nx3=nh;nx4=naz
    else
    	error("style not recognized.")
    end

    stream = open(join([list[1] ".seish"]))
    seek(stream, header_count["d1"])
    dt = read(stream,Float32)
    close(stream)
    d = zeros(Float32,nt,1)
    h = Array(Header,1)
    h[1] = InitSeisHeader()

    j = 1    
    for ix1 = 1 : nx1
    for ix2 = 1 : nx2
    for ix3 = 1 : nx3
    for ix4 = 1 : nx4
		h[1].tracenum = convert(typeof(h[1].tracenum),j)
		h[1].o1 = convert(typeof(h[1].o1),0)
		h[1].n1 = convert(typeof(h[1].n1),nt)
		h[1].d1 = convert(typeof(h[1].d1),dt)
		if (style=="sxsygxgy")
            h[1].isx = convert(typeof(h[1].isx),ix1 - 1 + min_isx)
			h[1].isy = convert(typeof(h[1].isy),ix2 - 1 + min_isy)
			h[1].igx = convert(typeof(h[1].igx),ix3 - 1 + min_igx)
			h[1].igy = convert(typeof(h[1].igy),ix4 - 1 + min_igy)
    		sx_rot = convert(Float32,(ix1 - 1 + min_isx)*dsx + osx);
    		sy_rot = convert(Float32,(ix2 - 1 + min_isy)*dsy + osy);
    		gx_rot = convert(Float32,(ix3 - 1 + min_igx)*dgx + ogx);
    		gy_rot = convert(Float32,(ix4 - 1 + min_igy)*dgy + ogy);
    		h[1].sx =  (sx_rot-osx)*cos(ang2) + (sy_rot-osy)*sin(ang2) + osx;
    		h[1].sy = -(sx_rot-osx)*sin(ang2) + (sy_rot-osy)*cos(ang2) + osy;
    		h[1].gx =  (gx_rot-ogx)*cos(ang2) + (gy_rot-ogy)*sin(ang2) + ogx;
    		h[1].gy = -(gx_rot-ogx)*sin(ang2) + (gy_rot-ogy)*cos(ang2) + ogy;		
        	h[1].hx = h[1].gx - h[1].sx
        	h[1].hy = h[1].gy - h[1].sy
			h[1].h = sqrt((h[1].hx^2) + (h[1].hy^2))
        	h[1].az = rad2deg*atan2((h[1].gy-h[1].sy),(h[1].gx-h[1].sx))
        	if (h[1].az < 0) 
	    		h[1].az += 360.0
	    	end
        	h[1].mx = h[1].sx + h[1].hx/(1 + gammainv);
        	h[1].my = h[1].sy + h[1].hy/(1 + gammainv);
        	mx_rot = (h[1].mx-omx)*cos(ang2) - (h[1].my-omy)*sin(ang2) + omx;
        	my_rot = (h[1].mx-omx)*sin(ang2) + (h[1].my-omy)*cos(ang2) + omy;
        	hx_rot = (h[1].hx-ohx)*cos(ang2) - (h[1].hy-ohy)*sin(ang2) + ohx;
        	hy_rot = (h[1].hx-ohx)*sin(ang2) + (h[1].hy-ohy)*cos(ang2) + ohy;
        	h[1].imx = convert(Int32,round((mx_rot-omx)/dmx))
        	h[1].imy = convert(Int32,round((my_rot-omy)/dmy))
        	h[1].ihx = convert(Int32,round((hx_rot-ohx)/dhx))
        	h[1].ihy = convert(Int32,round((hy_rot-ohy)/dhy))
        	h[1].ih = convert(Int32,round((h[1].h-oh)/dh))
        	h[1].iaz = convert(Int32,round((h[1].az-oaz)/daz))
		elseif (style=="mxmyhxhy")
            h[1].imx = convert(typeof(h[1].imx),ix1 - 1 + min_imx)
			h[1].imy = convert(typeof(h[1].imy),ix2 - 1 + min_imy)
			h[1].ihx = convert(typeof(h[1].ihx),ix3 - 1 + min_ihx)
			h[1].ihy = convert(typeof(h[1].ihy),ix4 - 1 + min_ihy)
    		mx_rot = convert(Float32,(ix1 - 1 + min_imx)*dmx + omx);
    		my_rot = convert(Float32,(ix2 - 1 + min_imy)*dmy + omy);
    		hx_rot = convert(Float32,(ix3 - 1 + min_ihx)*dhx + ohx);
    		hy_rot = convert(Float32,(ix4 - 1 + min_ihy)*dhy + ohy);
    		h[1].mx =  (mx_rot-omx)*cos(ang2) + (my_rot-omy)*sin(ang2) + omx;
    		h[1].my = -(mx_rot-omx)*sin(ang2) + (my_rot-omy)*cos(ang2) + omy;
    		h[1].hx =  (hx_rot-ohx)*cos(ang2) + (hy_rot-ohy)*sin(ang2) + ohx;
    		h[1].hy = -(hx_rot-ohx)*sin(ang2) + (hy_rot-ohy)*cos(ang2) + ohy;
            h[1].sx = h[1].mx - h[1].hx/(1 + gammainv);
            h[1].sy = h[1].my - h[1].hy/(1 + gammainv);
            h[1].gx = h[1].mx + h[1].hx*(1-(1/(1 + gammainv)));
            h[1].gy = h[1].my + h[1].hy*(1-(1/(1 + gammainv)));
    		sx_rot = (h[1].sx-osx)*cos(ang2) - (h[1].sy-osy)*sin(ang2) + osx;
    		sy_rot = (h[1].sx-osx)*sin(ang2) + (h[1].sy-osy)*cos(ang2) + osy;
    		gx_rot = (h[1].gx-ogx)*cos(ang2) - (h[1].gy-ogy)*sin(ang2) + ogx;
    		gy_rot = (h[1].gx-ogx)*sin(ang2) + (h[1].gy-ogy)*cos(ang2) + ogy;
            h[1].isx = convert(Int32,round((sx_rot-osx)/dsx))
            h[1].isy = convert(Int32,round((sy_rot-osy)/dsy))
            h[1].igx = convert(Int32,round((gx_rot-ogx)/dgx))
            h[1].igy = convert(Int32,round((gy_rot-ogy)/dgy))
            h[1].h = sqrt((h[1].hx^2) + (h[1].hy^2))
            h[1].az = rad2deg*atan2((h[1].gy-h[1].sy),(h[1].gx-h[1].sx))
            if (h[1].az < 0)
            	h[1].az += 360.0
            end
        	h[1].ih = convert(Int32,round((h[1].h-oh)/dh))
        	h[1].iaz = convert(Int32,round((h[1].az-oaz)/daz))
		elseif (style=="mxmyhaz")
			error("mxmyhaz not developed yet.")
		elseif (style=="sxsyhxhy")
			error("sxsyhxhy not developed yet.")
		elseif (style=="gxgyhxhy")
			error("gxgyhxhy not developed yet.")
		elseif (style=="sxsyhaz")
			error("sxsyhaz not developed yet.")
		elseif (style=="gxgyhaz")
			error("gxgyhaz not developed yet.")
		end
		h[1].selev = convert(typeof(h[1].selev),0)
		h[1].gelev = convert(typeof(h[1].gelev),0)
		h[1].trid = convert(typeof(h[1].trid),0)
    	SeisWrite(out,d,h,j)
    	j += 1
    end
    end
    end
    end

    out_d = join([out ".seisd"])
    out_h = join([out ".seish"])
    stream_d = open(out_d,"a+")
    stream_h = open(out_h,"a+")
	for ipatch = 1 : length(list)
		d_patch,h_patch,status = SeisRead(list[ipatch])
		#println("size(d_patch)=",size(d_patch))
		nx_patch = size(d_patch,2)
		nt_patch = h_patch[1].n1
		ot_patch = h_patch[1].o1
    	        min_it_patch = int(floor(ot_patch/dt)+1)
    	        max_it_patch = int(floor(ot_patch/dt)+nt_patch)
		min_ix1_patch = getfield(h_patch[1],symbol(key[2]))
		max_ix1_patch = getfield(h_patch[nx_patch],symbol(key[2]))
		min_ix2_patch = getfield(h_patch[1],symbol(key[3]))
		max_ix2_patch = getfield(h_patch[nx_patch],symbol(key[3]))
		min_ix3_patch = getfield(h_patch[1],symbol(key[4]))
		max_ix3_patch = getfield(h_patch[nx_patch],symbol(key[4]))
		min_ix4_patch = getfield(h_patch[1],symbol(key[5]))
		max_ix4_patch = getfield(h_patch[nx_patch],symbol(key[5]))
		
    		min_it_patch  >       1 ? tapti  =  it_WO : tapti  = 0
    		max_it_patch  <      nt ? taptf  =  it_WO : taptf  = 0    
    		min_ix1_patch > min_ix1 ? tapx1i = ix1_WO : tapx1i = 0
    		max_ix1_patch < max_ix1 ? tapx1f = ix1_WO : tapx1f = 0
    		min_ix2_patch > min_ix2 ? tapx2i = ix2_WO : tapx2i = 0
    		max_ix2_patch < max_ix2 ? tapx2f = ix2_WO : tapx2f = 0
    		min_ix3_patch > min_ix3 ? tapx3i = ix3_WO : tapx3i = 0
    		max_ix3_patch < max_ix3 ? tapx3f = ix3_WO : tapx3f = 0
    		min_ix4_patch > min_ix4 ? tapx4i = ix4_WO : tapx4i = 0
    		max_ix4_patch < max_ix4 ? tapx4f = ix4_WO : tapx4f = 0


	#println("min_it_patch=",min_it_patch)
	#println("max_it_patch=",max_it_patch)

		d_patch = Taper(d_patch,
				max_it_patch  -  min_it_patch + 1,
				max_ix1_patch - min_ix1_patch + 1,
				max_ix2_patch - min_ix2_patch + 1,
				max_ix3_patch - min_ix3_patch + 1,
				max_ix4_patch - min_ix4_patch + 1,
				tapti,taptf,tapx1i,tapx1f,tapx2i,tapx2f,tapx3i,tapx3f,tapx4i,tapx4f)
		
    	for j = 1 : nx_patch
    		h = h_patch[j]
        	if (style=="sxsygxgy")
    			itrace = (h.isx - min_isx)*nx2*nx3*nx4 + (h.isy - min_isy)*nx3*nx4 + (h.igx - min_igx)*nx4 + h.igy - min_igy + 1
			elseif (style=="mxmyhxhy")
    			itrace = (h.imx - min_imx)*nx2*nx3*nx4 + (h.imy - min_imy)*nx3*nx4 + (h.ihx - min_ihx)*nx4 + h.ihy - min_ihy + 1
			elseif (style=="mxmyhaz")
    			itrace = (h.imx - min_imx)*nx2*nx3*nx4 + (h.imy - min_imy)*nx3*nx4 + (h.ih - min_ih)*nx4 + h.iaz - min_iaz + 1
			elseif (style=="sxsyhxhy")
    			itrace = (h.isx - min_isx)*nx2*nx3*nx4 + (h.isy - min_isy)*nx3*nx4 + (h.ihx - min_ihx)*nx4 + h.ihy - min_ihy + 1
			elseif (style=="gxgyhxhy")
    			itrace = (h.igx - min_igx)*nx2*nx3*nx4 + (h.igy - min_igy)*nx3*nx4 + (h.ihx - min_ihx)*nx4 + h.ihy - min_ihy + 1
			elseif (style=="sxsyhaz")
    			itrace = (h.isx - min_isx)*nx2*nx3*nx4 + (h.isy - min_isy)*nx3*nx4 + (h.ih - min_ih)*nx4 + h.iaz - min_iaz + 1
			elseif (style=="gxgyhaz")
    			itrace = (h.igx - min_igx)*nx2*nx3*nx4 + (h.igy - min_igy)*nx3*nx4 + (h.ih - min_ih)*nx4 + h.iaz - min_iaz + 1
			end
			h.tracenum = convert(Int32,itrace)
    		position_d = 4*nt*(itrace - 1)
    		seek(stream_d,position_d)
    		d = read(stream_d,Float32,nt)
    		d[min_it_patch:max_it_patch] += d_patch[:,j]
    		seek(stream_d,position_d)
			write(stream_d,convert(Array{Float32,1},d[:]))
    		PutHeader(stream_h,h,itrace)
    	end
	end
    close(stream_d)
    close(stream_h)

end
