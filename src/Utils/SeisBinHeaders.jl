"""
    SeisBinHeaders(in,out; <keyword arguments>)

Sequentially bin seismic headers using the available grid information.

Keyword arguments should be consistent with SeisGeometry keyword arguments.

# Arguments
* `in`: filename of input, irregularly sampled data
* `out`: filename of output, regularly sampled data

# Keyword arguments
* `style="sxsygxgy"`: bin style. Options: "mxmyhxhy","mxmyhaz","sxsyhxhy","gxgyhxhy","sxsyhaz","gxgyhaz"
* `ang=90`: inline direction measured in degrees CC from East
* `gamma=1`: vp/vs ratio for PS Asymptotic Conversion Point gathers (use gamma=1 for PP data)
* `osx=0`,`osy=0`,`ogx=0`,`ogy=0` : origin for source and receiver coordinate system
* `omx=0`,`omy=0`,`ohx=0`,`ohy=0`: origin for midpoint and offset coordinate system
* `oaz=0`,`oh=0` : origin for azimuth and offset coordinate system
* `dsx=1`,`dsy=1`,`dgx=1`,`dgy=1`: source and receiver step-size
* `dmx=1`,`dmy=1`,`dhx=1`,`dhy=1`: midpoint and offset step-size
* `dh=1`,`daz=1`: offset and azimuth step-size
* `min_isx=0`,`max_isx=0`,`min_isy=0`,`max_isy=0`: grid extreme values for sources
* `min_igx=0`,`max_igx=0`,`min_igy=0`,`max_igy=0`: grid extreme values for receivers
* `min_imx=0`,`max_imx=0`,`min_imy=0`,`max_imy=0`: grid extreme values for midpoints
* `min_ihx=0`,`max_ihx=0`,`min_ihy=0`,`max_ihy=0`: grid extreme values for offsets
* `min_ih=0`,`max_ih=0`,`min_iaz=0`,`max_iaz=0`: grid extreme values for azimuth and offset
* `ntrace=10000`: maximum number of traces processed at a time

# Output
In file `out`, binned headers are created.

*Credits: Aaron Stanton,2017*

"""

function SeisBinHeaders(in,out;style="sxsygxgy",ang=90,gamma=1,osx=0,osy=0,ogx=0,ogy=0,omx=0,omy=0,ohx=0,ohy=0,oh=0,oaz=0,dsx=1,dsy=1,dgx=1,dgy=1,dmx=1,dmy=1,dhx=1,dhy=1,dh=1,daz=1,min_isx=0,max_isx=0,min_isy=0,max_isy=0,min_igx=0,max_igx=0,min_igy=0,max_igy=0,min_imx=0,max_imx=0,min_imy=0,max_imy=0,min_ihx=0,max_ihx=0,min_ihy=0,max_ihy=0,min_ih=0,max_ih=0,min_iaz=0,max_iaz=0,ntrace=10000)



rad2deg = 180/pi;
deg2rad = pi/180;
gammainv = 1/gamma;
if (ang > 90)
	ang2=-deg2rad*(ang-90)
else
	ang2=deg2rad*(90-ang)
end

naz=convert(Int32,360/daz)

if (style=="sxsygxgy")
	nsx = max_isx - min_isx + 1
	nsy = max_isy - min_isy + 1
	ngx = max_igx - min_igx + 1
	ngy = max_igy - min_igy + 1
	nx1=nsx;nx2=nsy;nx3=ngx;nx4=ngy;
	ox1=osx;ox2=osy;ox3=ogx;ox4=ogy;
	dx1=dsx;dx2=dsy;dx3=dgx;dx4=dgy;
	label2="sx";label3="sy";label4="gx";label5="gy";
	unit2="m";unit3="m";unit4="m";unit5="m";
elseif (style=="mxmyhxhy")
	nmx = max_imx - min_imx + 1
	nmy = max_imy - min_imy + 1
	nhx = max_ihx - min_ihx + 1
	nhy = max_ihy - min_ihy + 1
	nx1=nmx;nx2=nmy;nx3=nhx;nx4=nhy;
	ox1=omx;ox2=omy;ox3=ohx;ox4=ohy;
	dx1=dmx;dx2=dmy;dx3=dhx;dx4=dhy;
	label2="mx";label3="my";label4="hx";label5="hy";
	unit2="m";unit3="m";unit4="m";unit5="m";
elseif (style=="mxmyhaz")
	nmx = max_imx - min_imx + 1
	nmy = max_imy - min_imy + 1
	nh  = max_ih  - min_ih  + 1
	naz = max_iaz - min_iaz + 1
	nx1=nmx;nx2=nmy;nx3=nh;nx4=naz;
	ox1=omx;ox2=omy;ox3=oh;ox4=oaz;
	dx1=dmx;dx2=dmy;dx3=dh;dx4=daz;
	label2="mx";label3="my";label4="h";label5="az";
	unit2="m";unit3="m";unit4="m";unit5="Degrees";
elseif (style=="sxsyhxhy")
	nsx = max_isx - min_isx + 1
	nsy = max_isy - min_isy + 1
	nhx = max_ihx - min_ihx + 1
	nhy = max_ihy - min_ihy + 1
	nx1=nsx;nx2=nsy;nx3=nhx;nx4=nhy;
	ox1=osx;ox2=osy;ox3=ohx;ox4=ohy;
	dx1=dsx;dx2=dsy;dx3=dhx;dx4=dhy;
	label2="sx";label3="sy";label4="hx";label5="hy";
	unit2="m";unit3="m";unit4="m";unit5="m";
elseif (style=="gxgyhxhy")
	ngx = max_igx - min_igx + 1
	ngy = max_igy - min_igy + 1
	nhx = max_ihx - min_ihx + 1
	nhy = max_ihy - min_ihy + 1
	nx1=ngx;nx2=ngy;nx3=nhx;nx4=nhy;
	ox1=ogx;ox2=ogy;ox3=ohx;ox4=ohy;
	dx1=dgx;dx2=dgy;dx3=dhx;dx4=dhy;
	label2="gx";label3="gy";label4="hx";label5="hy";
	unit2="m";unit3="m";unit4="m";unit5="m";

elseif (style=="sxsyhaz")
	nsx = max_isx - min_isx + 1
	nsy = max_isy - min_isy + 1
	nh  = max_ih  - min_ih  + 1
	naz = max_iaz - min_iaz + 1
	nx1=nsx;nx2=nsy;nx3=nh;nx4=naz;
	ox1=osx;ox2=osy;ox3=oh;ox4=oaz;
	dx1=dsx;dx2=dsy;dx3=dh;dx4=daz;
	label2="sx";label3="sy";label4="hx";label5="az";
	unit2="m";unit3="m";unit4="m";unit5="Degrees";

elseif (style=="gxgyhaz")
 	ngx = max_igx - min_igx + 1
	ngy = max_igy - min_igy + 1
	nh  = max_ih  - min_ih  + 1
	naz = max_iaz - min_iaz + 1
	nx1=ngx;nx2=ngy;nx3=nh;nx4=naz;
	ox1=ogx;ox2=ogy;ox3=oh;ox4=oaz;
	dx1=dgx;dx2=dgy;dx3=dh;dx4=daz;
	label2="gx";label3="gy";label4="h";label5="az";
	unit2="m";unit3="m";unit4="m";unit5="Degrees";
else
	error("style not recognized.")
end



nx_out = nx1*nx2*nx3*nx4
stream_in = open(ParseHeaderName(in))

@compat nx_in = round(Int,filesize(stream_in)/(4*length(fieldnames(Header))))

seek(stream_in, header_count["n1"])
nt = read(stream_in,Int32)
println("nt= ",nt)
println("The final binned cube will have an approximate size of ", round(nx1*nx2*nx3*nx4*nt*4*1e-9,3)," Gb")
seek(stream_in, header_count["o1"])
o1 = read(stream_in,Float32)
seek(stream_in, header_count["d1"])
dt = read(stream_in,Float32)
close(stream_in)

extent = Extent(convert(Int32,nt),convert(Int32,nx1),convert(Int32,nx2),convert(Int32,nx3),convert(Int32,nx4),
	   convert(Float32,o1),convert(Float32,ox1),convert(Float32,ox2),convert(Float32,ox3),convert(Float32,ox4),
	   convert(Float32,dt),convert(Float32,dx1),convert(Float32,dx2),convert(Float32,dx3),convert(Float32,dx4),
	   "Time",label5,label4,label3,label2,
	   "s",unit5,unit4,unit3,unit2,
	   "")

DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
filename_d = join([DATAPATH out "@data@"])
filename_h = join([DATAPATH out "@headers@"])

#writes extent file
WriteTextHeader(out,extent,"native_float", 4,filename_d,filename_h)


h = Array{Header}(1)
h[1] = InitSeisHeader()
stream_out = open(filename_h,"a+")

j = 1

if (style=="sxsygxgy")

   for ix4 = 1 : nx4
       for ix3 = 1 : nx3
       	  for ix2 = 1 : nx2
		for ix1 = 1 : nx1
	       	  	h[1].tracenum = convert(typeof(h[1].tracenum),j)
			h[1].o1 = convert(typeof(h[1].o1),o1)
			h[1].n1 = convert(typeof(h[1].n1),nt)
			h[1].d1 = convert(typeof(h[1].d1),dt)
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
			h[1].iaz = convert(Int32,round((h[1].az-oaz)/daz))<naz?convert(Int32,round((h[1].az-oaz)/daz)):0
			h[1].selev = convert(typeof(h[1].selev),0)
			h[1].gelev = convert(typeof(h[1].gelev),0)
			h[1].trid = convert(typeof(h[1].trid),0)
			PutHeader(stream_out,h[1],j)
			j += 1
		end
     	  end
      end
   end

   h_in=SeisReadHeaders(in)
   hout = SeisReadHeaders(out)

   for k = 1 : nx_in

	itrace = (h_in[k].igy - min_igy)*nx3*nx2*nx1 + (h_in[k].igx - min_igx)*nx2*nx1 + (h_in[k].isy - min_isy)*nx1 + h_in[k].isx - min_isx + 1

	if (itrace > 0 && itrace <= nx_out)
	    	if (h_in[k].isx >= min_isx && h_in[k].isx <= max_isx && h_in[k].isy >= min_isy && h_in[k].isy <= max_isy && h_in[k].igx >= min_igx && h_in[k].igx <= max_igx && h_in[k].igy >= min_igy && h_in[k].igy <= max_igy)

			hout=GrabHeader(stream_out,itrace)

		     	hout.sx=h_in[k].sx
			hout.sy=h_in[k].sy
			hout.gx=h_in[k].gx
			hout.gy=h_in[k].gy
			hout.mx=h_in[k].mx
			hout.my=h_in[k].my
			hout.hx=h_in[k].hx
			hout.hy=h_in[k].hy
			hout.h=h_in[k].h
			hout.az=h_in[k].az
			hout.ang=h_in[k].ang

			PutHeader(stream_out,hout,itrace)
		end
	end
    end

 elseif (style=="mxmyhxhy")

    for ix4 = 1 : nx4
       for ix3 = 1 : nx3
       	  for ix2 = 1 : nx2
		for ix1 = 1 : nx1

			h[1].tracenum = convert(typeof(h[1].tracenum),j)
			h[1].o1 = convert(typeof(h[1].o1),o1)
			h[1].n1 = convert(typeof(h[1].n1),nt)
			h[1].d1 = convert(typeof(h[1].d1),dt)
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
			h[1].iaz = convert(Int32,round((h[1].az-oaz)/daz))<naz?convert(Int32,round((h[1].az-oaz)/daz)):0
			h[1].selev = convert(typeof(h[1].selev),0)
			h[1].gelev = convert(typeof(h[1].gelev),0)
			h[1].trid = convert(typeof(h[1].trid),0)
			PutHeader(stream_out,h[1],j)

			j += 1
		 end
     	 end
     end
  end

  h_in=SeisReadHeaders(in)
  hout = SeisReadHeaders(out)

  for k = 1 : nx_in

	itrace = (h_in[k].ihy - min_ihy)*nx3*nx2*nx1 + (h_in[k].ihx - min_ihx)*nx2*nx1 + (h_in[k].imy - min_imy)*nx1 + h_in[k].imx - min_imx + 1
	if (itrace > 0 && itrace <= nx_out)
	 if (h_in[k].imx >= min_imx && h_in[k].imx <= max_imx && h_in[k].imy >= min_imy && h_in[k].imy <= max_imy && h_in[k].ihx >= min_ihx && h_in[k].ihx <= max_ihx && h_in[k].ihy >= min_ihy && h_in[k].ihy <= max_ihy)

		 hout=GrabHeader(stream_out,itrace)
		 hout.sx=h_in[k].sx
		 hout.sy=h_in[k].sy
		 hout.gx=h_in[k].gx
		 hout.gy=h_in[k].gy
		 hout.mx=h_in[k].mx
		 hout.my=h_in[k].my
		 hout.hx=h_in[k].hx
		 hout.hy=h_in[k].hy
		 hout.h=h_in[k].h
		 hout.az=h_in[k].az
		 hout.ang=h_in[k].ang
		 PutHeader(stream_out,hout,itrace)

	 end
       end
 end

elseif (style=="mxmyhaz")
	for ix4 = 1 : nx4
         for ix3 = 1 : nx3
	  for ix2 = 1 : nx2
	   for ix1 = 1 : nx1

       	 	h[1].tracenum = convert(typeof(h[1].tracenum),j)
 		h[1].o1 = convert(typeof(h[1].o1),o1)
		h[1].n1 = convert(typeof(h[1].n1),nt)
		h[1].d1 = convert(typeof(h[1].d1),dt)
		h[1].imx = convert(typeof(h[1].imx),ix1 - 1 + min_imx)
		h[1].imy = convert(typeof(h[1].imy),ix2 - 1 + min_imy)
		h[1].ih = convert(typeof(h[1].ih),  ix3 - 1 + min_ih)
		h[1].iaz = convert(typeof(h[1].iaz),ix4 - 1 + min_iaz)
		mx_rot = convert(Float32,(ix1 - 1 + min_imx)*dmx + omx);
		my_rot = convert(Float32,(ix2 - 1 + min_imy)*dmy + omy);
		h[1].mx =  (mx_rot-omx)*cos(ang2) + (my_rot-omy)*sin(ang2) + omx;
		h[1].my = -(mx_rot-omx)*sin(ang2) + (my_rot-omy)*cos(ang2) + omy;
		h[1].h = convert(Float32,(ix3 - 1 + min_ih)*dh + oh);
		h[1].az = rad2deg*atan2((h[1].gy-h[1].sy),(h[1].gx-h[1].sx))
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
		h[1].selev = convert(typeof(h[1].selev),0)
		h[1].gelev = convert(typeof(h[1].gelev),0)
		h[1].trid = convert(typeof(h[1].trid),0)
		PutHeader(stream_out,h[1],j)

		j += 1
	end
     end
    end
   end

 h_in=SeisReadHeaders(in)
 hout = SeisReadHeaders(out)


 for k = 1 : nx_in
	itrace = (h_in[k].iaz - min_iaz)*nx3*nx2*nx1 + (h_in[k].ih - min_ih)*nx2*nx1 + (h_in[k].imy - min_imy)*nx1 + h_in[k].imx - min_imx + 1
	 if (itrace > 0 && itrace <= nx_out)
		 if (h_in[k].imx >= min_imx && h_in[k].imx <= max_imx && h_in[k].imy >= min_imy && h_in[k].imy <= max_imy && h_in[k].ih >= min_ih && h_in[k].ih <= max_ih && h_in[k].iaz >= min_iaz && h_in[k].iaz <= max_iaz)

			 hout[itrace].sx=h_in[k].sx
			 hout[itrace].sy=h_in[k].sy
			 hout[itrace].gx=h_in[k].gx
			 hout[itrace].gy=h_in[k].gy
			 hout[itrace].mx=h_in[k].mx
			 hout[itrace].my=h_in[k].my
			 hout[itrace].hx=h_in[k].hx
			 hout[itrace].hy=h_in[k].hy
			 hout[itrace].h=h_in[k].h
			 hout[itrace].az=h_in[k].az
			 hout[itrace].ang=h_in[k].ang
			 PutHeader(stream_out,hout[itrace],itrace)

		 end
	 end
 end


elseif (style=="sxsyhxhy")
 for ix4 = 1 : nx4
  for ix3 = 1 : nx3
   for ix2 = 1 : nx2
    for ix1 = 1 : nx1

	h[1].tracenum = convert(typeof(h[1].tracenum),j)
	h[1].o1 = convert(typeof(h[1].o1),o1)
	h[1].n1 = convert(typeof(h[1].n1),nt)
	h[1].d1 = convert(typeof(h[1].d1),dt)
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
	h[1].iaz = convert(Int32,round((h[1].az-oaz)/daz))<naz?convert(Int32,round((h[1].az-oaz)/daz)):0
	h[1].selev = convert(typeof(h[1].selev),0)
	h[1].gelev = convert(typeof(h[1].gelev),0)
	h[1].trid = convert(typeof(h[1].trid),0)
	PutHeader(stream_out,h[1],j)

	j += 1
    end
   end
  end
 end


 h_in=SeisReadHeaders(in)
 hout = SeisReadHeaders(out)

 for k = 1 : nx_in
  itrace = (h_in[k].ihy - min_ihy)*nx3*nx2*nx1 + (h_in[k].ihx - min_ihx)*nx2*nx1 + (h_in[k].isy - min_isy)*nx1 + h_in[k].isx - min_isx + 1
	if (itrace > 0 && itrace <= nx_out)

	 if (h_in[k].isx >= min_isx && h_in[k].isx <= max_isx && h_in[k].isy >= min_isy && h_in[k].isy <= max_isy && h_in[k].ihx >= min_ihx && h_in[k].ihx <= max_ihx && h_in[k].ihy >= min_ihy && h_in[k].ihy <= max_ihy)

	 	hout[itrace].sx=h_in[k].sx
		hout[itrace].sy=h_in[k].sy
		hout[itrace].gx=h_in[k].gx
		hout[itrace].gy=h_in[k].gy
		hout[itrace].mx=h_in[k].mx
		hout[itrace].my=h_in[k].my
		hout[itrace].hx=h_in[k].hx
		hout[itrace].hy=h_in[k].hy
		hout[itrace].h=h_in[k].h
		hout[itrace].az=h_in[k].az
		hout[itrace].ang=h_in[k].ang
                PutHeader(stream_out,hout[itrace],itrace)

         end
      end
  end



elseif (style=="gxgyhxhy")
 for ix4 = 1 : nx4
  for ix3 = 1 : nx3
   for ix2 = 1 : nx2
    for ix1 = 1 : nx1

	h[1].tracenum = convert(typeof(h[1].tracenum),j)
	h[1].o1 = convert(typeof(h[1].o1),o1)
	h[1].n1 = convert(typeof(h[1].n1),nt)
	h[1].d1 = convert(typeof(h[1].d1),dt)
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
	h[1].iaz = convert(Int32,round((h[1].az-oaz)/daz))<naz?convert(Int32,round((h[1].az-oaz)/daz)):0
	h[1].selev = convert(typeof(h[1].selev),0)
	h[1].gelev = convert(typeof(h[1].gelev),0)
	h[1].trid = convert(typeof(h[1].trid),0)
	PutHeader(stream_out,h[1],j)

	j += 1
     end
    end
   end
  end

 h_in=SeisReadHeaders(in)
 hout = SeisReadHeaders(out)

 for k = 1 : nx_in
	itrace = (h_in[k].ihy - min_ihy)*nx3*nx2*nx1 + (h_in[k].ihx - min_ihx)*nx2*nx1 + (h_in[k].igy - min_igy)*nx1 + h_in[k].igx - min_igx + 1
	if (itrace > 0 && itrace <= nx_out)

	 if (h_in[k].igx >= min_igx && h_in[k].igx <= max_igx && h_in[k].igy >= min_igy && h_in[k].igy <= max_igy && h_in[k].ihx >= min_ihx && h_in[k].ihx <= max_ihx && h_in[k].ihy >= min_ihy && h_in[k].ihy <= max_ihy)

		 hout[itrace].sx=h_in[k].sx
		 hout[itrace].sy=h_in[k].sy
		 hout[itrace].gx=h_in[k].gx
		 hout[itrace].gy=h_in[k].gy
		 hout[itrace].mx=h_in[k].mx
		 hout[itrace].my=h_in[k].my
		 hout[itrace].hx=h_in[k].hx
		 hout[itrace].hy=h_in[k].hy
		 hout[itrace].h=h_in[k].h
		 hout[itrace].az=h_in[k].az
		 hout[itrace].ang=h_in[k].ang
		 PutHeader(stream_out,hout[itrace],itrace)

	 end
       end
     end

elseif (style=="sxsyhaz")
	error("sxsyhaz not developed yet.")

elseif (style=="gxgyhaz")
	error("gxgyhaz not developed yet.")

end

close(stream_out)

end
