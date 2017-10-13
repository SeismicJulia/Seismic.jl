"""
    SeisBinData(in,out; <keyword arguments>)

Sequentially bin seismic data using already binned trace headers (SeisBinHeaders).
Input arguments should be consistent with SeisBinHeaders input arguments.

# Arguments
* `in::AbstractString`: filename of input, irregularly sampled data
* `out::AbstractString`: filename of output, regularly sampled data

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
In file `out`, the binned data is created.

*Credits: Aaron Stanton, 2017*

"""

function SeisBinData(in::AbstractString,out::AbstractString;style="sxsygxgy",ang=90,gamma=1,osx=0,osy=0,ogx=0,ogy=0,omx=0,omy=0,ohx=0,ohy=0,oh=0,oaz=0,dsx=1,dsy=1,dgx=1,dgy=1,dmx=1,dmy=1,dhx=1,dhy=1,dh=1,daz=1,min_isx=0,max_isx=0,min_isy=0,max_isy=0,min_igx=0,max_igx=0,min_igy=0,max_igy=0,min_imx=0,max_imx=0,min_imy=0,max_imy=0,min_ihx=0,max_ihx=0,min_ihy=0,max_ihy=0,min_ih=0,max_ih=0,min_iaz=0,max_iaz=0,ntrace=10000)


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

	elseif (style=="mxmyhxhy")
		nmx = max_imx - min_imx + 1
		nmy = max_imy - min_imy + 1
		nhx = max_ihx - min_ihx + 1
		nhy = max_ihy - min_ihy + 1

		nx1=nmx;nx2=nmy;nx3=nhx;nx4=nhy;

	elseif (style=="mxmyhaz")
	  nmx = max_imx - min_imx + 1
		nmy = max_imy - min_imy + 1
		nh  = max_ih  - min_ih  + 1
		naz = max_iaz - min_iaz + 1

		nx1=nmx;nx2=nmy;nx3=nh;nx4=naz;

	elseif (style=="sxsyhxhy")
	  nsx = max_isx - min_isx + 1
	  nsy = max_isy - min_isy + 1
		nhx = max_ihx - min_ihx + 1
		nhy = max_ihy - min_ihy + 1

		nx1=nsx;nx2=nsy;nx3=nhx;nx4=nhy;

	elseif (style=="gxgyhxhy")
		ngx = max_igx - min_igx + 1
		ngy = max_igy - min_igy + 1
		nhx = max_ihx - min_ihx + 1
		nhy = max_ihy - min_ihy + 1

		nx1=ngx;nx2=ngy;nx3=nhx;nx4=nhy;

	elseif (style=="sxsyhaz")
	  nsx = max_isx - min_isx + 1
	 	nsy = max_isy - min_isy + 1
		nh  = max_ih  - min_ih  + 1
		naz = max_iaz - min_iaz + 1

		nx1=nsx;nx2=nsy;nx3=nh;nx4=naz;

	elseif (style=="gxgyhaz")
    ngx = max_igx - min_igx + 1
		ngy = max_igy - min_igy + 1
		nh  = max_ih  - min_ih  + 1
		naz = max_iaz - min_iaz + 1

		nx1=ngx;nx2=ngy;nx3=nh;nx4=naz;

	else
		error("style not recognized.")
	end

	nx_out = nx1*nx2*nx3*nx4

	stream = open(ParseHeaderName(in))

	@compat nx_in = round(Int,filesize(stream)/(4*length(fieldnames(Header))))

	seek(stream, header_count["n1"])
	nt = read(stream,Int32)

	seek(stream, header_count["o1"])
	o1 = read(stream,Float32)

	seek(stream, header_count["d1"])
	dt = read(stream,Float32)

	close(stream)

	d = zeros(Float32,nt,1)

	out_d = ParseDataName(out)
	out_h = ParseHeaderName(out)
	stream_d = open(out_d,"a")


#Create 0 volume
	j = 1
	for ix4 = 1 : nx4
			for ix3 = 1 : nx3
				 for ix2 = 1 : nx2
					 for ix1 = 1 : nx1

		    		position_d = 4*nt*(j - 1)
		    		seek(stream_d,position_d)
						write(stream_d,d[:,1])
						j += 1
		    	end
				end
  	  end
	end
	close(stream_d)

	stream_d = open(out_d,"a")

	if (style=="sxsygxgy")
		j = 1
		while j <= nx_in
			d,h,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = ntrace
			for k = 1 : num_traces_in
				itrace = (h[k].igy - min_igy)*nx1*nx2*nx3 + (h[k].igx - min_igx)*nx2*nx1 + (h[k].isy - min_isy)*nx1 + h[k].igx - min_igx + 1
				if (itrace > 0 && itrace <= nx_out)

					if (h[k].isx >= min_isx && h[k].isx <= max_isx && h[k].isy >= min_isy && h[k].isy <= max_isy && h[k].igx >= min_igx && h[k].igx <= max_igx && h[k].igy >= min_igy && h[k].igy <= max_igy)

						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])

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
				itrace = (h[k].ihy - min_ihy)*nx3*nx2*nx1 + (h[k].ihx - min_ihx)*nx2*nx1 + (h[k].imy - min_imy)*nx1 + h[k].imx - min_imx + 1
	 		 if (itrace > 0 && itrace <= nx_out)

	 				 if (h[k].imx >= min_imx && h[k].imx <= max_imx && h[k].imy >= min_imy && h[k].imy <= max_imy && h[k].ihx >= min_ihx && h[k].ihx <= max_ihx && h[k].ihy >= min_ihy && h[k].ihy <= max_ihy)


						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])

					end
				end
			end
			j += num_traces_in
		end

	elseif (style=="mxmyhaz")
		j = 1
		while j <= nx_in
			d,h_in,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)

			for k = 1 : num_traces_in

				itrace = (h_in[k].iaz - min_iaz)*nx3*nx2*nx1 + (h_in[k].ih - min_ih)*nx2*nx1 + (h_in[k].imy - min_imy)*nx1 + h_in[k].imx - min_imx + 1
			  if (itrace > 0 && itrace <= nx_out)

			 		 if (h_in[k].imx >= min_imx && h_in[k].imx <= max_imx && h_in[k].imy >= min_imy && h_in[k].imy <= max_imy && h_in[k].ih >= min_ih && h_in[k].ih <= max_ih && h_in[k].iaz >= min_iaz && h_in[k].iaz <= max_iaz)

						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])

					end
				end
			end
			j += num_traces_in
		end

	elseif (style=="sxsyhxhy")
		j = 1
		while j <= nx_in
			d,h_in,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)

			for k = 1 : num_traces_in
				itrace = (h_in[k].ihy - min_ihy)*nx3*nx2*nx1 + (h_in[k].ihx - min_ihx)*nx2*nx1 + (h_in[k].isy - min_isy)*nx1 + h_in[k].isx - min_isx + 1
				if (itrace > 0 && itrace <= nx_out)

					 if (h_in[k].isx >= min_isx && h_in[k].isx <= max_isx && h_in[k].isy >= min_isy && h_in[k].isy <= max_isy && h_in[k].ihx >= min_ihx && h_in[k].ihx <= max_ihx && h_in[k].ihy >= min_ihy && h_in[k].ihy <= max_ihy)

						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])
					end
				end
			end
			j += num_traces_in
		end

	elseif (style=="gxgyhxhy")
		j = 1
		while j <= nx_in
			d,h_in,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)
			for k = 1 : num_traces_in
				itrace = (h_in[k].ihy - min_ihy)*nx3*nx2*nx1 + (h_in[k].ihx - min_ihx)*nx2*nx1 + (h_in[k].igy - min_igy)*nx1 + h_in[k].igx - min_igx + 1
				if (itrace > 0 && itrace <= nx_out)

					 if (h_in[k].igx >= min_igx && h_in[k].igx <= max_igx && h_in[k].igy >= min_igy && h_in[k].igy <= max_igy && h_in[k].ihx >= min_ihx && h_in[k].ihx <= max_ihx && h_in[k].ihy >= min_ihy && h_in[k].ihy <= max_ihy)

						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])

					end
				end
			end
			j += num_traces_in
		end


	elseif (style=="sxsyhaz")
		j = 1
		while j <= nx_in
			d,h_in,e = SeisRead(in,group="some",itrace=j,ntrace=ntrace)
			num_traces_in = size(d[:,:],2)
			for k = 1 : num_traces_in
				itrace = (h_in[k].iaz - min_iaz)*nx3*nx2*nx1 + (h_in[k].ih - min_ih)*nx2*nx1 + (h_in[k].isy - min_isy)*nx1 + h_in[k].isx - min_isx + 1
				if (itrace > 0 && itrace <= nx_out)

					 if (h_in[k].isx >= min_isx && h_in[k].isx <= max_isx && h_in[k].isy >= min_isy && h_in[k].isy <= max_isy && h_in[k].ih >= min_ih && h_in[k].ih <= max_ih && h_in[k].ias >= min_iaz && h_in[k].iaz <= max_iaz)

						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])
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

				itrace = (h[k].iaz - min_iaz)*nx1*nx2*nx3 + (h[k].ih - min_ih)*nx1*nx2 + (h[k].igy - min_igy)*nx1 + h[k].igx - min_igx + 1
				if (itrace > 0 && itrace <= nx_out)

					if (h[k].igx >= min_igx && h[k].igx <= max_igx && h[k].igy >= min_igy && h[k].igy <= max_igy && h[k].ih >= min_ih && h[k].ih <= max_ih && h[k].iaz >= min_iaz && h[k].iaz <= max_iaz)

						position_d = 4*nt*(itrace - 1)
						seek(stream_d,position_d)
						write(stream_d,d[:,k])

					end
				end
			end
			j += num_traces_in
		end
	end
	close(stream_d)

end
