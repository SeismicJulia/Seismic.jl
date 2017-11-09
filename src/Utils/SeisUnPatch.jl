"""
    SeisUnPatch(in,out;<keyword arguments>)

Reconstruct a 5D data volume from a set of 5D data patches.

# Arguments
* `in::Array{AbstractString,1}`: array containing filename of patches
* `out::AbstractString`: filename for reconstructed volume

# Keyword arguments
* `style="sxsygxgy"`: bin style. Options: "mxmyhxhy","mxmyhaz","sxsyhxhy","gxgyhxhy","sxsyhaz","gxgyhaz"
* `min_isx=0`,`max_isx=0`,`min_isy=0`,`max_isy=0`: grid extreme values for sources
* `min_igx=0`,`max_igx=0`,`min_igy=0`,`max_igy=0`: grid extreme values for receivers
* `min_imx=0`,`max_imx=0`,`min_imy=0`,`max_imy=0`: grid extreme values for midpoints
* `min_ihx=0`,`max_ihx=0`,`min_ihy=0`,`max_ihy=0`: grid extreme values for offsets
* `min_ih=0`,`max_ih=0`,`min_iaz=0`,`max_iaz=0`: grid extreme values for azimuth and offset
* `it_WL=9e9`,`it_WO=0` : length and overlapping samples in time patches
* `ix1_WL=9e9`,`ix1_WO=0`:length and overlapping samples in first space dimension
* `ix2_WL=9e9`,`ix2_WO=0`,`ix3_WL=9e9`,`ix3_WO=0`,`ix4_WL=9e9`,`ix4_WO=0`
* `nt=0`: time samples of reconstructed cube
* `ang=90`: inline direction measured in degrees CC from East
* `gamma=1`: vp/vs ratio for PS Asymptotic Conversion Point gathers (use gamma=1 for PP data)
* `osx=0`,`osy=0`,`ogx=0`,`ogy=0` : origin for source and receiver coordinate system
* `omx=0`,`omy=0`,`ohx=0`,`ohy=0`: origin for midpoint and offset coordinate system
* `oaz=0`,`oh=0` : origin for azimuth and offset coordinate system
* `dsx=1`,`dsy=1`,`dgx=1`,`dgy=1`: source and receiver step-size
* `dmx=1`,`dmy=1`,`dhx=1`,`dhy=1`: midpoint and offset step-size
* `dh=1`,`daz=1`: offset and azimuth step-size

# Output
In file `out`, the 5D reconstructed volume is created.

*Credits: A. Stanton, F Carozzi, 2017*
"""

function SeisUnPatch(patch_names::Array{AbstractString,1},out::AbstractString;style="sxsygxgy",min_isx=0,max_isx=0,min_isy=0,max_isy=0,min_igx=0,max_igx=0,min_igy=0,max_igy=0,min_imx=0,max_imx=0,min_imy=0,max_imy=0,min_ihx=0,max_ihx=0,min_ihy=0,max_ihy=0,min_ih=0,max_ih=0,min_iaz=0,max_iaz=0,it_WL=9e9,it_WO=0,ix1_WL=9e9,ix1_WO=0,ix2_WL=9e9,ix2_WO=0,ix3_WL=9e9,ix3_WO=0,ix4_WL=9e9,ix4_WO=0,nt=0,ang=90,gamma=1,osx=0,osy=0,ogx=0,ogy=0,omx=0,omy=0,ohx=0,ohy=0,oh=0,oaz=0,dsx=1,dsy=1,dgx=1,dgy=1,dmx=1,dmy=1,dhx=1,dhy=1,dh=1,daz=1)


    if (style == "sxsygxgy")
	key = ["t","isx","isy","igx","igy"]
	min_ix1 = min_isx
	max_ix1 = max_isx
	min_ix2 = min_isy
	max_ix2 = max_isy
	min_ix3 = min_igx
	max_ix3 = max_igx
	min_ix4 = min_igy
	max_ix4 = max_igy
    elseif (style=="mxmyhxhy")
	key = ["t","imx","imy","ihx","ihy"]
	min_ix1 = min_imx
	max_ix1 = max_imx
	min_ix2 = min_imy
	max_ix2 = max_imy
	min_ix3 = min_ihx
	max_ix3 = max_ihx
	min_ix4 = min_ihy
	max_ix4 = max_ihy
    elseif (style=="mxmyhaz")
	key = ["t","imx","imy","ih","iaz"]
	min_ix1 = min_imx
	max_ix1 = max_imx
	min_ix2 = min_imy
	max_ix2 = max_imy
	min_ix3 = min_ih
	max_ix3 = max_ih
	min_ix4 = min_iaz
	max_ix4 = max_iaz
    elseif (style=="sxsyhxhy")
	key = ["t","isx","isy","ihx","ihy"]
	min_ix1 = min_isx
	max_ix1 = max_isx
	min_ix2 = min_isy
	max_ix2 = max_isy
	min_ix3 = min_ihx
	max_ix3 = max_ihx
	min_ix4 = min_ihy
	max_ix4 = max_ihy
    elseif (style=="gxgyhxhy")
	key = ["t","igx","igy","ihx","ihy"]
	min_ix1 = min_igx
	max_ix1 = max_igx
	min_ix2 = min_igy
	max_ix2 = max_igy
	min_ix3 = min_ihx
	max_ix3 = max_ihx
	min_ix4 = min_ihy
	max_ix4 = max_ihy
    elseif (style=="sxsyhaz")
	key = ["t","isx","isy","ih","iaz"]
	min_ix1 = min_isx
	max_ix1 = max_isx
	min_ix2 = min_isy
	max_ix2 = max_isy
	min_ix3 = min_ih
	max_ix3 = max_ih
	min_ix4 = min_iaz
	max_ix4 = max_iaz
    elseif (style=="gxgyhaz")
	key = ["t","igx","igy","ih","iaz"]
	min_ix1 = min_igx
	max_ix1 = max_igx
	min_ix2 = min_igy
		max_ix2 = max_igy
	min_ix3 = min_ih
	max_ix3 = max_ih
	min_ix4 = min_iaz
	max_ix4 = max_iaz
    else
	error("style not defined.")
    end
    nx1 = max_ix1 - min_ix1 + 1
    nx2 = max_ix2 - min_ix2 + 1
    nx3 = max_ix3 - min_ix3 + 1
    nx4 = max_ix4 - min_ix4 + 1

    it_WL  = it_WL  > nt  ? nt  : it_WL
    ix1_WL = ix1_WL > nx1 ? nx1 : ix1_WL
	ix2_WL = ix2_WL > nx2 ? nx2 : ix2_WL
    ix3_WL = ix3_WL > nx3 ? nx3 : ix3_WL
    ix4_WL = ix4_WL > nx4 ? nx4 : ix4_WL

    rad2deg = 180/pi
    deg2rad = pi/180
    gammainv = 1/gamma
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
    nh = max_ih - min_ih + 1
    naz = max_iaz - min_iaz + 1

    if (style=="sxsygxgy")
	nx1=nsx; nx2=nsy; nx3=ngx; nx4=ngy
    elseif (style=="mxmyhxhy")
	nx1=nmx; nx2=nmy; nx3=nhx; nx4=nhy
    elseif (style=="mxmyhaz")
	nx1=nmx; nx2=nmy; nx3=nh; nx4=naz
    elseif (style=="sxsyhxhy")
	nx1=nsx; nx2=nsy; nx3=nhx; nx4=nhy
    elseif (style=="gxgyhxhy")
	nx1=ngx; nx2=ngy; nx3=nhx; nx4=nhy
    elseif (style=="sxsyhaz")
	nx1=nsx; nx2=nsy; nx3=nh; nx4=naz
    elseif (style=="gxgyhaz")
	nx1=ngx; nx2=ngy; nx3=nh; nx4=naz
	else
	    error("style not recognized.")
	end


    filename_data = ParseDataName(patch_names[1])
    filename_headers = ParseHeaderName(patch_names[1])
    extent = ReadTextHeader(patch_names[1])


    dt = extent.d1
    ot = extent.o1
    d = zeros(Float32,nt,1)
    h = Array{Header}(1)
    h[1] = InitSeisHeader()

    extent = Seismic.Extent(nt, max_ix1-min_ix1+1, max_ix2-min_ix2+1,
                            max_ix3-min_ix3+1, max_ix4-min_ix4+1,convert(Float32,ot), 0, 0, 0,
                            0, convert(Float32,dt), 1, 1, 1, 1, "", "", "", "", "", "", "", "",
                            "", "", "")
println("creating headers")
   j = 1

     for ix1 = 1 : nx1
       for ix2 = 1 : nx2
         for ix3 = 1 : nx3
           for ix4 = 1 : nx4

              h[1].tracenum = convert(typeof(h[1].tracenum),j)

        h[1].o1 = convert(typeof(h[1].o1),ot)
		    h[1].n1 = convert(typeof(h[1].n1),nt)
		    h[1].d1 = convert(typeof(h[1].d1),dt)

		    if (style=="sxsygxgy")
			h[1].isx = convert(typeof(h[1].isx),ix1 - 1 + min_isx)
			h[1].isy = convert(typeof(h[1].isy),ix2 - 1 + min_isy)
			h[1].igx = convert(typeof(h[1].igx),ix3 - 1 + min_igx)
			h[1].igy = convert(typeof(h[1].igy),ix4 - 1 + min_igy)
			sx_rot = convert(Float32,(ix1 - 1 + min_isx)*dsx + osx)
			sy_rot = convert(Float32,(ix2 - 1 + min_isy)*dsy + osy)
			gx_rot = convert(Float32,(ix3 - 1 + min_igx)*dgx + ogx)
			gy_rot = convert(Float32,(ix4 - 1 + min_igy)*dgy + ogy)
			h[1].sx =  (sx_rot-osx)*cos(ang2) + (sy_rot-osy)*sin(ang2) + osx
			h[1].sy = -(sx_rot-osx)*sin(ang2) + (sy_rot-osy)*cos(ang2) + osy
			h[1].gx =  (gx_rot-ogx)*cos(ang2) + (gy_rot-ogy)*sin(ang2) + ogx
			h[1].gy = -(gx_rot-ogx)*sin(ang2) + (gy_rot-ogy)*cos(ang2) + ogy
			h[1].hx = h[1].gx - h[1].sx
			h[1].hy = h[1].gy - h[1].sy
			h[1].h = sqrt((h[1].hx^2) + (h[1].hy^2))
			h[1].az = rad2deg*atan2((h[1].gy-h[1].sy),(h[1].gx-h[1].sx))
			if (h[1].az < 0)
			    h[1].az += 360.0
			end
			h[1].mx = h[1].sx + h[1].hx/(1 + gammainv)
			h[1].my = h[1].sy + h[1].hy/(1 + gammainv)
			mx_rot = (h[1].mx-omx)*cos(ang2) - (h[1].my-omy)*sin(ang2) + omx
			my_rot = (h[1].mx-omx)*sin(ang2) + (h[1].my-omy)*cos(ang2) + omy
			hx_rot = (h[1].hx-ohx)*cos(ang2) - (h[1].hy-ohy)*sin(ang2) + ohx
			hy_rot = (h[1].hx-ohx)*sin(ang2) + (h[1].hy-ohy)*cos(ang2) + ohy
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
			mx_rot = convert(Float32,(ix1 - 1 + min_imx)*dmx + omx)
			my_rot = convert(Float32,(ix2 - 1 + min_imy)*dmy + omy)
			hx_rot = convert(Float32,(ix3 - 1 + min_ihx)*dhx + ohx)
			hy_rot = convert(Float32,(ix4 - 1 + min_ihy)*dhy + ohy)
			h[1].mx =  (mx_rot-omx)*cos(ang2) + (my_rot-omy)*sin(ang2) + omx
			h[1].my = -(mx_rot-omx)*sin(ang2) + (my_rot-omy)*cos(ang2) + omy
			h[1].hx =  (hx_rot-ohx)*cos(ang2) + (hy_rot-ohy)*sin(ang2) + ohx
			h[1].hy = -(hx_rot-ohx)*sin(ang2) + (hy_rot-ohy)*cos(ang2) + ohy
			h[1].sx = h[1].mx - h[1].hx/(1 + gammainv)
			h[1].sy = h[1].my - h[1].hy/(1 + gammainv)
			h[1].gx = h[1].mx + h[1].hx*(1-(1/(1 + gammainv)))
			h[1].gy = h[1].my + h[1].hy*(1-(1/(1 + gammainv)))
			sx_rot = (h[1].sx-osx)*cos(ang2) - (h[1].sy-osy)*sin(ang2) + osx
			sy_rot = (h[1].sx-osx)*sin(ang2) + (h[1].sy-osy)*cos(ang2) + osy
			gx_rot = (h[1].gx-ogx)*cos(ang2) - (h[1].gy-ogy)*sin(ang2) + ogx
			gy_rot = (h[1].gx-ogx)*sin(ang2) + (h[1].gy-ogy)*cos(ang2) + ogy
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
			mx_rot = convert(Float32,(ix1 - 1 + min_imx)*dmx + omx)
			my_rot = convert(Float32,(ix2 - 1 + min_imy)*dmy + omy)
			h[1].mx =  (mx_rot-omx)*cos(ang2) + (my_rot-omy)*sin(ang2) + omx
			h[1].my = -(mx_rot-omx)*sin(ang2) + (my_rot-omy)*cos(ang2) + omy
			h[1].h = convert(Float32,(ix3 - 1 + min_ih)*dh + oh)
			h[1].az = convert(Float32,(ix4 - 1 + min_iaz)*daz + oaz)
			if (h[1].az <= 90)
			    h[1].hx = h[1].h*cos(deg2rad*h[1].az)
			    h[1].hy = h[1].h*sin(deg2rad*h[1].az)
			elseif (h[1].az > 90 && h[1].az <= 180)
			    h[1].hx =-h[1].h*cos(pi-(deg2rad*h[1].az))
			    h[1].hy = h[1].h*sin(pi-(deg2rad*h[1].az))
			elseif (h[1].az > 180 && h[1].az <= 270)
			    h[1].hx =-h[1].h*cos((deg2rad*h[1].az)-pi)
			    h[1].hy =-h[1].h*sin((deg2rad*h[1].az)-pi)
			else
			    h[1].hx = h[1].h*cos(2*pi-(deg2rad*h[1].az))
			    h[1].hy =-h[1].h*sin(2*pi-(deg2rad*h[1].az))
			end
			h[1].sx = h[1].mx - h[1].hx/(1 + gammainv)
			h[1].sy = h[1].my - h[1].hy/(1 + gammainv)
			h[1].gx = h[1].mx + h[1].hx*(1-(1/(1 + gammainv)))
			h[1].gy = h[1].my + h[1].hy*(1-(1/(1 + gammainv)))
			sx_rot = (h[1].sx-osx)*cos(ang2) - (h[1].sy-osy)*sin(ang2) + osx
			sy_rot = (h[1].sx-osx)*sin(ang2) + (h[1].sy-osy)*cos(ang2) + osy
			gx_rot = (h[1].gx-ogx)*cos(ang2) - (h[1].gy-ogy)*sin(ang2) + ogx
			gy_rot = (h[1].gx-ogx)*sin(ang2) + (h[1].gy-ogy)*cos(ang2) + ogy
			h[1].isx = convert(Int32,round((sx_rot-osx)/dsx))
			h[1].isy = convert(Int32,round((sy_rot-osy)/dsy))
			h[1].igx = convert(Int32,round((gx_rot-ogx)/dgx))
			h[1].igy = convert(Int32,round((gy_rot-ogy)/dgy))
			hx_rot = (h[1].hx-ohx)*cos(ang2) - (h[1].hy-ohy)*sin(ang2) + ohx
			hy_rot = (h[1].hx-ohx)*sin(ang2) + (h[1].hy-ohy)*cos(ang2) + ohy
			h[1].ihx = convert(Int32,round((hx_rot-ohx)/dhx))
			h[1].ihy = convert(Int32,round((hy_rot-ohy)/dhy))
		    elseif (style=="sxsyhxhy")
			h[1].isx = convert(typeof(h[1].isx),ix1 - 1 + min_isx)
			h[1].isy = convert(typeof(h[1].isy),ix2 - 1 + min_isy)
			h[1].ihx = convert(typeof(h[1].ihx),ix3 - 1 + min_ihx)
			h[1].ihy = convert(typeof(h[1].ihy),ix4 - 1 + min_ihy)
			sx_rot = convert(Float32,(ix1 - 1 + min_isx)*dsx + osx)
			sy_rot = convert(Float32,(ix2 - 1 + min_isy)*dsy + osy)
			hx_rot = convert(Float32,(ix3 - 1 + min_ihx)*dhx + ohx)
			hy_rot = convert(Float32,(ix4 - 1 + min_ihy)*dhy + ohy)
			h[1].sx =  (sx_rot-osx)*cos(ang2) + (sy_rot-osy)*sin(ang2) + osx
			h[1].sy = -(sx_rot-osx)*sin(ang2) + (sy_rot-osy)*cos(ang2) + osy
			h[1].hx =  (hx_rot-ohx)*cos(ang2) + (hy_rot-ohy)*sin(ang2) + ohx
			h[1].hy = -(hx_rot-ohx)*sin(ang2) + (hy_rot-ohy)*cos(ang2) + ohy
			h[1].gx = h[1].sx + h[1].hx
			h[1].gy = h[1].sy + h[1].hy
			h[1].h = sqrt((h[1].hx^2) + (h[1].hy^2))
			h[1].az = rad2deg*atan2((h[1].gy-h[1].sy),(h[1].gx-h[1].sx))
			if (h[1].az < 0)
			    h[1].az += 360.0
			end
			h[1].mx = h[1].sx + h[1].hx/(1 + gammainv)
			h[1].my = h[1].sy + h[1].hy/(1 + gammainv)
			mx_rot = (h[1].mx-omx)*cos(ang2) - (h[1].my-omy)*sin(ang2) + omx
			my_rot = (h[1].mx-omx)*sin(ang2) + (h[1].my-omy)*cos(ang2) + omy
			hx_rot = (h[1].hx-ohx)*cos(ang2) - (h[1].hy-ohy)*sin(ang2) + ohx
			hy_rot = (h[1].hx-ohx)*sin(ang2) + (h[1].hy-ohy)*cos(ang2) + ohy
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
			gx_rot = convert(Float32,(ix1 - 1 + min_igx)*dgx + ogx)
			gy_rot = convert(Float32,(ix2 - 1 + min_igy)*dgy + ogy)
			hx_rot = convert(Float32,(ix3 - 1 + min_ihx)*dhx + ohx)
			hy_rot = convert(Float32,(ix4 - 1 + min_ihy)*dhy + ohy)
			h[1].gx =  (gx_rot-ogx)*cos(ang2) + (gy_rot-ogy)*sin(ang2) + ogx
			h[1].gy = -(gx_rot-ogx)*sin(ang2) + (gy_rot-ogy)*cos(ang2) + ogy
			h[1].hx =  (hx_rot-ohx)*cos(ang2) + (hy_rot-ohy)*sin(ang2) + ohx
			h[1].hy = -(hx_rot-ohx)*sin(ang2) + (hy_rot-ohy)*cos(ang2) + ohy
			h[1].sx = h[1].gx - h[1].hx
			h[1].sy = h[1].gy - h[1].hy
			h[1].h = sqrt((h[1].hx^2) + (h[1].hy^2))
			h[1].az = rad2deg*atan2((h[1].gy-h[1].sy),(h[1].gx-h[1].sx))
			if (h[1].az < 0)
			    h[1].az += 360.0
			end
			h[1].mx = h[1].sx + h[1].hx/(1 + gammainv)
			h[1].my = h[1].sy + h[1].hy/(1 + gammainv)
			mx_rot = (h[1].mx-omx)*cos(ang2) - (h[1].my-omy)*sin(ang2) + omx
			my_rot = (h[1].mx-omx)*sin(ang2) + (h[1].my-omy)*cos(ang2) + omy
			hx_rot = (h[1].hx-ohx)*cos(ang2) - (h[1].hy-ohy)*sin(ang2) + ohx
			hy_rot = (h[1].hx-ohx)*sin(ang2) + (h[1].hy-ohy)*cos(ang2) + ohy
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


    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
    out_d = join([DATAPATH out "@data@"])
    out_h = join([DATAPATH out "@headers@"])


    for ipatch = 1 : length(patch_names)

   println("ipatch= ",ipatch," ",patch_names[ipatch])
	stream_d = open(out_d,"a+")
	stream_h = open(out_h,"a+")

	d_patch,h_patch,e_patch = SeisRead(patch_names[ipatch])

	nx_patch=e_patch.n2*e_patch.n3*e_patch.n4*e_patch.n5

	nt_patch = h_patch[1].n1

	d_patch=reshape(d_patch,Int(nt_patch),Int(nx_patch))

	ot_patch = h_patch[1].o1
	min_it_patch = Int(floor((ot_patch-ot)/dt)+1)
	max_it_patch = Int(floor((ot_patch-ot)/dt)+nt_patch)


	min_ix1_patch = getfield(h_patch[1],Symbol(key[2]))
	max_ix1_patch = getfield(h_patch[nx_patch],Symbol(key[2]))
	min_ix2_patch = getfield(h_patch[1],Symbol(key[3]))
	max_ix2_patch = getfield(h_patch[nx_patch],Symbol(key[3]))
	min_ix3_patch = getfield(h_patch[1],Symbol(key[4]))
	max_ix3_patch = getfield(h_patch[nx_patch],Symbol(key[4]))
	min_ix4_patch = getfield(h_patch[1],Symbol(key[5]))
	max_ix4_patch = getfield(h_patch[nx_patch],Symbol(key[5]))


	tapti  =  min_it_patch >       1 ? it_WO  : 0
	taptf  =  max_it_patch <      nt ? it_WO  : 0
	tapx1i = min_ix1_patch > min_ix1 ? ix1_WO : 0
	tapx1f = max_ix1_patch < max_ix1 ? ix1_WO : 0
	tapx2i = min_ix2_patch > min_ix2 ? ix2_WO : 0
	tapx2f = max_ix2_patch < max_ix2 ? ix2_WO : 0
	tapx3i = min_ix3_patch > min_ix3 ? ix3_WO : 0
	tapx3f = max_ix3_patch < max_ix3 ? ix3_WO : 0
	tapx4i = min_ix4_patch > min_ix4 ? ix4_WO : 0
	tapx4f = max_ix4_patch < max_ix4 ? ix4_WO : 0


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
        itrace = (h.igy - min_igy)*nx3*nx2*nx1 + (h.igx - min_igx)*nx2*nx1 + (h.isy - min_isy)*nx1 + h.isx - min_isx + 1
	    elseif (style=="mxmyhxhy")
          itrace = (h.ihy - min_ihy)*nx3*nx2*nx1 + (h.ihx - min_ihx)*nx2*nx1 + (h.imy - min_imy)*nx1 + h.imx - min_imx + 1
	    elseif (style=="mxmyhaz")
		      # itrace = (h.imx - min_imx)*nx2*nx3*nx4 + (h.imy - min_imy)*nx3*nx4 + (h.ih - min_ih)*nx4 + h.iaz - min_iaz + 1
          itrace=(h.iaz-min_iaz)*nx3*nx2*nx1 + (h.ih - min_ih)*nx2*nx1 + (h.imy - min_imy)*nx1 + h.imx - min_imx+1
	    elseif (style=="sxsyhxhy")
          itrace =(h.ihy-min_ihy)*nx3*nx2*nx1 + (h.ihx - min_ihx)*nx2*nx1 + (h.isy - min_isy)*nx1 + h.isx - min_isx+1
	    elseif (style=="gxgyhxhy")
        itrace =(h.ihy-min_ihy)*nx3*nx2*nx1 + (h.ihx - min_ihx)*nx2*nx1 + (h.igy - min_igy)*nx1 + h.igx - min_igx+1
	    elseif (style=="sxsyhaz")
        itrace =(h.iaz-min_iaz)*nx3*nx2*nx1 + (h.ih - min_ih)*nx2*nx1 + (h.isy - min_isy)*nx1 + h.isx - min_isx+1
	    elseif (style=="gxgyhaz")
        itrace =(h.iaz-min_iaz)*nx3*nx2*nx1 + (h.ih - min_ih)*nx2*nx1 + (h.igy - min_igy)*nx1 + h.igx - min_igx+1
      end

      if j==1
        println(min_it_patch," ",max_it_patch," ",nt," ",nx_patch," ",ipatch," ",size(d_patch))
      end

      h.tracenum = itrace
	    h.n1 = nt
	    h.o1 = ot
	    position_d = 4*nt*(itrace - 1)
	    seek(stream_d,position_d)
	    d = read(stream_d,Float32,nt)
	    d[min_it_patch:max_it_patch] += d_patch[:,j]


	    seek(stream_d,position_d)
	    write(stream_d,d)
	    PutHeader(stream_h,h,itrace)
	end
	close(stream_d)
	close(stream_h)
    end

end

function Taper(d,nt,nx1,nx2,nx3,nx4,tapti,taptf,tapx1i,tapx1f,tapx2i,tapx2f,tapx3i,tapx3f,tapx4i,tapx4f)
    tx1=1; tx2=1; tx3=1; tx4=1


    for ix1 = 1 : nx1

		if (ix1>=1   && ix1<=tapx1i)
	        tx1 = 1.0/(tapx1i-1)*(ix1-1)

	     end

	     if (ix1>tapx1i && ix1<=nx1-tapx1f)
	        tx1 = 1
	     end
	     if (ix1>nx1-tapx1f && ix1<=nx1)
	        tx1 = 1-1.0/(tapx1f-1)*(ix1-1-nx1+tapx1f)

	      end

	     for ix2 = 1 : nx2
	        if (ix2>=1   && ix2<=tapx2i)
		          tx2 = 1.0/(tapx2i-1)*(ix2-1)

	        end

	        if (ix2>tapx2i && ix2<=nx2-tapx2f)
		          tx2 = 1
	        end
	        if (ix2>nx2-tapx2f && ix2<=nx2)
		       tx2 = 1-1.0/(tapx2f-1)*(ix2-1-nx2+tapx2f)

		end


	        for ix3 = 1 : nx3
		          if (ix3>=1   && ix3<=tapx3i)
		              tx3 = 1.0/(tapx3i-1)*(ix3-1)

		          end

		          if (ix3>tapx3i && ix3<=nx3-tapx3f)
		              tx3 = 1
		          end
		          if (ix3>nx3-tapx3f && ix3<=nx3)
		              tx3 = 1-1.0/(tapx3f-1)*(ix3-1-nx3+tapx3f)

             end


		            for ix4 = 1 : nx4
		                if (ix4>=1   && ix4<=tapx4i)
			                   tx4 = 1.0/(tapx4i-1)*(ix4-1)

		                end


		                if (ix4>tapx4i && ix4<=nx4-tapx4f)
			                   tx4 = 1
		                end
		                if (ix4>nx4-tapx4f && ix4<=nx4)
                        		tx4 = 1-1.0/(tapx4f-1)*(ix4-1-nx4+tapx4f)

			      end


		                ix = (ix4-1)*nx2*nx3*nx1 + (ix3-1)*nx1*nx2 + (ix2-1)*nx1 + ix1
		                for it = 1 : nt
			                   if (it>=1 && it<=tapti)

			                        tt = 1.0/(tapti-1)*(it-1)
			                   end

			                   if (it>tapti && it<=nt-taptf)
			                        tt = 1
			                   end
			                   if (it>nt-taptf && it<=nt)
			                        tt = 1-1.0/(taptf-1)*(it-1-nt+taptf)

			                   end

			      d[it,ix] = tt*tx1*tx2*tx3*tx4*d[it,ix]
		          end
		     end
	       end
	  end
      end
 return d
end
