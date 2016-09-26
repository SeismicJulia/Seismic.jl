"""
**SeisBin**

*Bin seismic data onto regular grid using trace headers*

**IN**

* in: filename of input, irregularly sampled data
* out: filename of output, regularly sampled data
* style="sxsygxgy"
* ang=90
* gamma=1
* osx=0
* osy=0
* ogx=0
* ogy=0
* omx=0
* omy=0
* ohx=0
* ohy=0
* oh=0
* oaz=0
* dsx=1
* dsy=1
* dgx=1
* dgy=1
* dmx=1
* dmy=1
* dhx=1
* dhy=1
* dh=1
* daz=1
* min_isx=0
* max_isx=0
* min_isy=0
* max_isy=0
* min_igx=0
* max_igx=0
* min_igy=0
* max_igy=0
* min_imx=0
* max_imx=0
* min_imy=0
* max_imy=0
* min_ihx=0
* max_ihx=0
* min_ihy=0
* max_ihy=0
* min_ih=0
* max_ih=0
* min_iaz=0
* max_iaz=0

**OUT**

*Credits: Aaron Stanton, 2015*

"""

function SeisBin(in,out;style="sxsygxgy",ang=90,gamma=1,osx=0,osy=0,ogx=0,ogy=0,omx=0,omy=0,ohx=0,ohy=0,oh=0,oaz=0,dsx=1,dsy=1,dgx=1,dgy=1,dmx=1,dmy=1,dhx=1,dhy=1,dh=1,daz=1,min_isx=0,max_isx=0,min_isy=0,max_isy=0,min_igx=0,max_igx=0,min_igy=0,max_igy=0,min_imx=0,max_imx=0,min_imy=0,max_imy=0,min_ihx=0,max_ihx=0,min_ihy=0,max_ihy=0,min_ih=0,max_ih=0,min_iaz=0,max_iaz=0,ntrace=10000)

	rad2deg = 180/pi;
	deg2rad = pi/180;
	gammainv = 1/gamma;
	if (ang > 90) 
		ang2=-deg2rad*(ang-90)
	else 
		ang2=deg2rad*(90-ang)
	end

	nsx = max_isx - min_isx + 1
	nsy = max_isy - min_isy + 1
	ngx = max_igx - min_igx + 1
	ngy = max_igy - min_igy + 1
	nmx = max_imx - min_imx + 1
	nmy = max_imy - min_imy + 1
	nhx = max_ihx - min_ihx + 1
	nhy = max_ihy - min_ihy + 1
	nh  = max_ih  - min_ih  + 1
	naz = max_iaz - min_iaz + 1

	if (style=="sxsygxgy")
		nx1=nsx;nx2=nsy;nx3=ngx;nx4=ngy;
		ox1=osx;ox2=osy;ox3=ogx;ox4=ogy;
		dx1=dsx;dx2=dsy;dx3=dgx;dx4=dgy;
		label2="sx";label3="sy";label4="gx";label5="gy";
		unit2="m";unit3="m";unit4="m";unit5="m";
	elseif (style=="mxmyhxhy")
		nx1=nmx;nx2=nmy;nx3=nhx;nx4=nhy;
		ox1=omx;ox2=omy;ox3=ohx;ox4=ohy;
		dx1=dmx;dx2=dmy;dx3=dhx;dx4=dhy;
		label2="mx";label3="my";label4="hx";label5="hy";
		unit2="m";unit3="m";unit4="m";unit5="m";
	elseif (style=="mxmyhaz")
		nx1=nmx;nx2=nmy;nx3=nh;nx4=naz;
		ox1=omx;ox2=omy;ox3=oh;ox4=oaz;
		dx1=dmx;dx2=dmy;dx3=dh;dx4=daz;
		label2="mx";label3="my";label4="h";label5="az";
		unit2="m";unit3="m";unit4="m";unit5="Degrees";
	elseif (style=="sxsyhxhy")
		nx1=nsx;nx2=nsy;nx3=nhx;nx4=nhy;
		ox1=osx;ox2=osy;ox3=ohx;ox4=ohy;
		dx1=dsx;dx2=dsy;dx3=dhx;dx4=dhy;
		label2="sx";label3="sy";label4="hx";label5="hy";
		unit2="m";unit3="m";unit4="m";unit5="m";
	elseif (style=="gxgyhxhy")
		nx1=ngx;nx2=ngy;nx3=nhx;nx4=nhy;
		ox1=ogx;ox2=ogy;ox3=ohx;ox4=ohy;
		dx1=dgx;dx2=dgy;dx3=dhx;dx4=dhy;
		label2="gx";label3="gy";label4="hx";label5="hy";
		unit2="m";unit3="m";unit4="m";unit5="m";
	elseif (style=="sxsyhaz")
		nx1=nsx;nx2=nsy;nx3=nh;nx4=naz;
		ox1=osx;ox2=osy;ox3=oh;ox4=oaz;
		dx1=dsx;dx2=dsy;dx3=dh;dx4=daz;
		label2="sx";label3="sy";label4="hx";label5="az";
		unit2="m";unit3="m";unit4="m";unit5="Degrees";
	elseif (style=="gxgyhaz")
		nx1=ngx;nx2=ngy;nx3=nh;nx4=naz;
		ox1=ogx;ox2=ogy;ox3=oh;ox4=oaz;
		dx1=dgx;dx2=dgy;dx3=dh;dx4=daz;
		label2="gx";label3="gy";label4="h";label5="az";
		unit2="m";unit3="m";unit4="m";unit5="Degrees";
	else
		error("style not recognized.")
	end
	nx_out = nx1*nx2*nx3*nx4

	stream = open(ParseHeaderName(in))
	@compat nx_in = int(filesize(stream)/(4*length(fieldnames(Header))))
	seek(stream, header_count["n1"])
	nt = read(stream,Int32)
	seek(stream, header_count["d1"])
	dt = read(stream,Float32)
	close(stream)
	d = zeros(Float32,nt,1)
	h = Array(Header,1)
	h[1] = InitSeisHeader()

	extent = Extent(convert(Int32,nt),convert(Int32,nx1),convert(Int32,nx2),convert(Int32,nx3),convert(Int32,nx4),
		   convert(Float32,0),convert(Float32,ox1),convert(Float32,ox2),convert(Float32,ox3),convert(Float32,ox4),
		   convert(Float32,dt),convert(Float32,dx1),convert(Float32,dx2),convert(Float32,dx3),convert(Float32,dx4),
		   "Time",label2,label3,label4,label5,
		   "s",unit2,unit3,unit4,unit5,
		   "")
	
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
						h[1].imx = convert(typeof(h[1].imx),ix1 - 1 + min_imx)
						h[1].imy = convert(typeof(h[1].imy),ix2 - 1 + min_imy)
						h[1].ih = convert(typeof(h[1].ih),  ix3 - 1 + min_ih)
						h[1].iaz = convert(typeof(h[1].iaz),ix4 - 1 + min_iaz)
						mx_rot = convert(Float32,(ix1 - 1 + min_imx)*dmx + omx);
						my_rot = convert(Float32,(ix2 - 1 + min_imy)*dmy + omy);
						h[1].mx =  (mx_rot-omx)*cos(ang2) + (my_rot-omy)*sin(ang2) + omx;
						h[1].my = -(mx_rot-omx)*sin(ang2) + (my_rot-omy)*cos(ang2) + omy;
						h[1].h = convert(Float32,(ix3 - 1 + min_ih)*dh + oh);
						h[1].az = convert(Float32,(ix4 - 1 + min_iaz)*daz + oaz);
						if (h[1].az <= 90)
							h[1].hx = h[1].h*cos(deg2rad*h[1].az);
							h[1].hy = h[1].h*sin(deg2rad*h[1].az);
						elseif (h[1].az > 90 && h[1].az <= 180)
							h[1].hx =-h[1].h*cos(pi-(deg2rad*h[1].az));
							h[1].hy = h[1].h*sin(pi-(deg2rad*h[1].az));
						elseif (h[1].az > 180 && h[1].az <= 270)
							h[1].hx =-h[1].h*cos((deg2rad*h[1].az)-pi);
							h[1].hy =-h[1].h*sin((deg2rad*h[1].az)-pi);
						else 
							h[1].hx = h[1].h*cos(2*pi-(deg2rad*h[1].az));
							h[1].hy =-h[1].h*sin(2*pi-(deg2rad*h[1].az));
						end					
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
						hx_rot = (h[1].hx-ohx)*cos(ang2) - (h[1].hy-ohy)*sin(ang2) + ohx;
						hy_rot = (h[1].hx-ohx)*sin(ang2) + (h[1].hy-ohy)*cos(ang2) + ohy;
						h[1].ihx = convert(Int32,round((hx_rot-ohx)/dhx))
						h[1].ihy = convert(Int32,round((hy_rot-ohy)/dhy))
					elseif (style=="sxsyhxhy")
						h[1].isx = convert(typeof(h[1].isx),ix1 - 1 + min_isx)
						h[1].isy = convert(typeof(h[1].isy),ix2 - 1 + min_isy)
						h[1].ihx = convert(typeof(h[1].ihx),ix3 - 1 + min_ihx)
						h[1].ihy = convert(typeof(h[1].ihy),ix4 - 1 + min_ihy)
						sx_rot = convert(Float32,(ix1 - 1 + min_isx)*dsx + osx);
						sy_rot = convert(Float32,(ix2 - 1 + min_isy)*dsy + osy);
						hx_rot = convert(Float32,(ix3 - 1 + min_ihx)*dhx + ohx);
						hy_rot = convert(Float32,(ix4 - 1 + min_ihy)*dhy + ohy);
						h[1].sx =  (sx_rot-osx)*cos(ang2) + (sy_rot-osy)*sin(ang2) + osx;
						h[1].sy = -(sx_rot-osx)*sin(ang2) + (sy_rot-osy)*cos(ang2) + osy;
						h[1].hx =  (hx_rot-ohx)*cos(ang2) + (hy_rot-ohy)*sin(ang2) + ohx;
						h[1].hy = -(hx_rot-ohx)*sin(ang2) + (hy_rot-ohy)*cos(ang2) + ohy;
						h[1].gx = h[1].sx + h[1].hx;
						h[1].gy = h[1].sy + h[1].hy;
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
					elseif (style=="gxgyhxhy")
						h[1].igx = convert(typeof(h[1].igx),ix1 - 1 + min_igx)
						h[1].igy = convert(typeof(h[1].igy),ix2 - 1 + min_igy)
						h[1].ihx = convert(typeof(h[1].ihx),ix3 - 1 + min_ihx)
						h[1].ihy = convert(typeof(h[1].ihy),ix4 - 1 + min_ihy)
						gx_rot = convert(Float32,(ix1 - 1 + min_igx)*dgx + ogx);
						gy_rot = convert(Float32,(ix2 - 1 + min_igy)*dgy + ogy);
						hx_rot = convert(Float32,(ix3 - 1 + min_ihx)*dhx + ohx);
						hy_rot = convert(Float32,(ix4 - 1 + min_ihy)*dhy + ohy);
						h[1].gx =  (gx_rot-ogx)*cos(ang2) + (gy_rot-ogy)*sin(ang2) + ogx;
						h[1].gy = -(gx_rot-ogx)*sin(ang2) + (gy_rot-ogy)*cos(ang2) + ogy;
						h[1].hx =  (hx_rot-ohx)*cos(ang2) + (hy_rot-ohy)*sin(ang2) + ohx;
						h[1].hy = -(hx_rot-ohx)*sin(ang2) + (hy_rot-ohy)*cos(ang2) + ohy;
						h[1].sx = h[1].gx - h[1].hx;
						h[1].sy = h[1].gy - h[1].hy;
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
					elseif (style=="sxsyhaz")
						error("sxsyhaz not developed yet.")
					elseif (style=="gxgyhaz")
						error("gxgyhaz not developed yet.")
					end
					h[1].selev = convert(typeof(h[1].selev),0)
					h[1].gelev = convert(typeof(h[1].gelev),0)
					h[1].trid = convert(typeof(h[1].trid),0)
					SeisWrite(out,d,h,extent,itrace=j)
					j += 1
				end
			end
		end
	end

	out_d = ParseDataName(out)
	out_h = ParseHeaderName(out)
	stream_d = open(out_d,"a")
	stream_h = open(out_h,"a")

	if (style=="sxsygxgy")
		j = 1
		while j <= nx_in	
			d,h,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)
			for k = 1 : num_traces_in
				itrace = (h[k].isx - min_isx)*nx2*nx3*nx4 + (h[k].isy - min_isy)*nx3*nx4 + (h[k].igx - min_igx)*nx4 + h[k].igy - min_igy + 1
				if (h[k].isx >= min_isx && h[k].isx <= max_isx && h[k].isy >= min_isy && h[k].isy <= max_isy && h[k].igx >= min_igx && h[k].igx <= max_igx && h[k].igy >= min_igy && h[k].igy <= max_igy)
					h[k].tracenum = convert(Int32,itrace)
					if (itrace > 0 && itrace <= nx_out)
						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])
						PutHeader(stream_h,h[k],itrace)
					end
				end
			end
			j += num_traces_in
		end
	elseif (style=="mxmyhxhy")
		j = 1
		while j <= nx_in	
			d,h,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)
			for k = 1 : num_traces_in
				itrace = (h[k].imx - min_imx)*nx2*nx3*nx4 + (h[k].imy - min_imy)*nx3*nx4 + (h[k].ihx - min_ihx)*nx4 + h[k].ihy - min_ihy + 1
				if (h[k].imx >= min_imx && h[k].imx <= max_imx && h[k].imy >= min_imy && h[k].imy <= max_imy && h[k].ihx >= min_ihx && h[k].ihx <= max_ihx && h[k].ihy >= min_ihy && h[k].ihy <= max_ihy)
					h[k].tracenum = convert(Int32,itrace)
					if (itrace > 0 && itrace <= nx_out)
						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])
						PutHeader(stream_h,h[k],itrace)
					end
				end
			end
			j += num_traces_in
		end
	elseif (style=="mxmyhaz")
		j = 1
		while j <= nx_in	
			d,h,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)
			for k = 1 : num_traces_in
				itrace = (h[k].imx - min_imx)*nx2*nx3*nx4 + (h[k].imy - min_imy)*nx3*nx4 + (h[k].ih - min_ih)*nx4 + h[k].iaz - min_iaz + 1
				if (h[k].imx >= min_imx && h[k].imx <= max_imx && h[k].imy >= min_imy && h[k].imy <= max_imy && h[k].ih >= min_ih && h[k].ih <= max_ih && h[k].iaz >= min_iaz && h[k].iaz <= max_iaz)
					h[k].tracenum = convert(Int32,itrace)
					if (itrace > 0 && itrace <= nx_out)
						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])
						PutHeader(stream_h,h[k],itrace)
					end
				end
			end
			j += num_traces_in
		end
	elseif (style=="sxsyhxhy")
		j = 1
		while j <= nx_in	
			d,h,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)
			for k = 1 : num_traces_in
				itrace = (h[k].isx - min_isx)*nx2*nx3*nx4 + (h[k].isy - min_isy)*nx3*nx4 + (h[k].ihx - min_ihx)*nx4 + h[k].ihy - min_ihy + 1
				if (h[k].isx >= min_isx && h[k].isx <= max_isx && h[k].isy >= min_isy && h[k].isy <= max_isy && h[k].ihx >= min_ihx && h[k].ihx <= max_ihx && h[k].ihy >= min_ihy && h[k].ihy <= max_ihy)
					h[k].tracenum = convert(Int32,itrace)
					if (itrace > 0 && itrace <= nx_out)
						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])
						PutHeader(stream_h,h[k],itrace)
					end
				end
			end
			j += num_traces_in
		end
	elseif (style=="gxgyhxhy")
		j = 1
		while j <= nx_in	
			d,h,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)
			for k = 1 : num_traces_in
				itrace = (h[k].igx - min_igx)*nx2*nx3*nx4 + (h[k].igy - min_igy)*nx3*nx4 + (h[k].ihx - min_ihx)*nx4 + h[k].ihy - min_ihy + 1
				if (h[k].igx >= min_igx && h[k].igx <= max_igx && h[k].igy >= min_igy && h[k].igy <= max_igy && h[k].ihx >= min_ihx && h[k].ihx <= max_ihx && h[k].ihy >= min_ihy && h[k].ihy <= max_ihy)
					h[k].tracenum = convert(Int32,itrace)
					if (itrace > 0 && itrace <= nx_out)
						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])
						PutHeader(stream_h,h[k],itrace)
					end
				end
			end
			j += num_traces_in
		end
	elseif (style=="sxsyhaz")
		j = 1
		while j <= nx_in	
			d,h,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)
			for k = 1 : num_traces_in
				itrace = (h[k].isx - min_isx)*nx2*nx3*nx4 + (h[k].isy - min_isy)*nx3*nx4 + (h[k].ih - min_ih)*nx4 + h[k].iaz - min_iaz + 1
				if (h[k].isx >= min_isx && h[k].isx <= max_isx && h[k].isy >= min_isy && h[k].isy <= max_isy && h[k].ih >= min_ih && h[k].ih <= max_ih && h[k].iaz >= min_iaz && h[k].iaz <= max_iaz)
					h[k].tracenum = convert(Int32,itrace)
					if (itrace > 0 && itrace <= nx_out)
						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])
						PutHeader(stream_h,h[k],itrace)
					end
				end
			end
			j += num_traces_in
		end
	elseif (style=="gxgyhaz")
		j = 1
		while j <= nx_in	
			d,h,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)
			for k = 1 : num_traces_in
				itrace = (h[k].igx - min_igx)*nx2*nx3*nx4 + (h[k].igy - min_igy)*nx3*nx4 + (h[k].ih - min_ih)*nx4 + h[k].iaz - min_iaz + 1
				if (h[k].igx >= min_igx && h[k].igx <= max_igx && h[k].igy >= min_igy && h[k].igy <= max_igy && h[k].ih >= min_ih && h[k].ih <= max_ih && h[k].iaz >= min_iaz && h[k].iaz <= max_iaz)
					h[k].tracenum = convert(Int32,itrace)
					if (itrace > 0 && itrace <= nx_out)
						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])
						PutHeader(stream_h,h[k],itrace)
					end
				end
			end
			j += num_traces_in
		end
	end	
	close(stream_d)
	close(stream_h)

end
