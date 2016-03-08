"""
**SeisPlot**

*Plot 2d seismic data with color, wiggles, or overlay. Optionally can print to a file. *

**IN**   
* d: 2d data 

* style="color"
* wiggle_fill_color="k"
* wiggle_line_color="k"
* wiggle_trace_increment=1
* xcur=1.2
* cmap="Greys"
* aspect="auto"
* pclip=98
* vmin="NULL"
* vmax="NULL"
* title=" "
* xlabel=" "
* xunits=" "
* ylabel=" "
* yunits=" "
* ox=0
* dx=1
* oy=0
* dy=1
* dpi=100
* wbox=6
* hbox=6
* name="NULL"
* interpolation="Hanning"
* titlesize=20
* labelsize=15
* ticksize=10
* fignum="NULL"
* plot_type="TX"
* fmax=100

**OUT**  

*Credits: Aaron Stanton, 2015*

"""
function SeisPlot(in;style="color",wiggle_fill_color="k",wiggle_line_color="k",wiggle_trace_increment=1,xcur=1.2,cmap="Greys",aspect="auto",pclip=98,vmin="NULL",vmax="NULL",title=" ",xlabel=" ",xunits=" ",ylabel=" ",yunits="",ox=0,dx=1,oy=0,dy=1,dpi=100,wbox=6,hbox=6,name="NULL",interpolation="Hanning",titlesize=20,labelsize=15,ticksize=10,fignum="NULL",plot_type="TX",fmax=100)

	if (vmin=="NULL" || vmax=="NULL")
		if (pclip<=100)
			a = -quantile(abs(in[:]),(pclip/100))
		else
			a = -quantile(abs(in[:]),1)*pclip/100
		end
		b = -a
	else
		a = vmin
		b = vmax
	end
	plt.ion()
	if (fignum == "NULL")
		fig = plt.figure(figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
	else
		fig = plt.figure(num=fignum,figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
	end
   if plot_type == "TX" # time space plot
	if (style != "wiggles")
		im = plt.imshow(in,cmap=cmap,vmin=a,vmax=b,extent=[ox,ox + (size(in,2)-1)*dx,oy + (size(in,1)-1)*dy,oy],aspect=aspect,interpolation=interpolation)
	end
	if (style != "color")
		y = oy+dy*collect(0:1:size(in,1)-1)
		x = ox+dx*collect(0:1:size(in,2)-1)
		delta = wiggle_trace_increment*dx
		hmin = minimum(x)
		hmax = maximum(x)
		dmax = maximum(in[:])
		alpha = xcur*delta/dmax
		for k = 1:wiggle_trace_increment:size(in,2)
			x_vert = Float64[]
			y_vert = Float64[]
			s = in[:,k]
			s[1]=0
			s[end]=0
			sp = (s+abs(s))/2
			im = plt.plot( s*alpha + x[k],y,wiggle_line_color)
			if (style != "overlay")
				plt.fill(sp*alpha + x[k],y,wiggle_fill_color)
			end
		end
		plt.axis([ox,ox + (size(in,2)-1)*dx,oy + (size(in,1)-1)*dy,oy])
	end
   elseif plot_type == "FK"
	xlabel = "Wavenumber"
	xunits = "(1/m)"
	ylabel = "Frequency"
	yunits = "Hz"
	dk = 1/dx/size(in[:,:],2)
	kmin = -dk*size(in[:,:],2)/2
	kmax =  dk*size(in[:,:],2)/2
	df = 1/dy/size(in[:,:],1)
	FMAX = df*size(in[:,:],1)/2 
	if fmax > FMAX
		fmax = FMAX
	end
	nf = convert(Int32,floor((size(in[:,:],1)/2)*fmax/FMAX))
	D = abs(fftshift(fft(in[:,:])))
	D = D[round(Int,end/2):round(Int,end/2)+nf,:]
	if (vmin=="NULL" || vmax=="NULL")
		a = 0.
		if (pclip<=100)
			b = quantile(abs(D[:]),(pclip/100))
		else
			b = quantile(abs(D[:]),1)*pclip/100
		end
	end
	im = plt.imshow(D,cmap=cmap,vmin=a,vmax=b,extent=[kmin,kmax,fmax,0],aspect=aspect,interpolation=interpolation)
   elseif plot_type == "Amplitude"
	xlabel = "Frequency"
	xunits = "(Hz)"
	ylabel = "Amplitude"
	yunits = ""
	nx = size(in[:,:],2)
	df = 1/dy/size(in[:,:],1)
	FMAX = df*size(in[:,:],1)/2 
	if fmax > FMAX
		fmax = FMAX
	end
	nf = convert(Int32,floor((size(in[:,:],1)/2)*fmax/FMAX))
	y = fftshift(sum(abs(fft(in[:,:],1)),2))/nx
	y = y[round(Int,end/2):round(Int,end/2)+nf]
	norm = maximum(y[:])
	if (norm > 0.)
		y = y/norm
	end
	x = collect(0:df:fmax)
	im = plt.plot(x,y)
	plt.title(title)
	plt.xlabel(join([xlabel " " xunits]))
	plt.ylabel(join([ylabel " " yunits]))
	plt.axis([0,fmax,0,1.1])
   else
	error("plot_type not recognized.")
   end 
	plt.title(title, fontsize=titlesize)
	plt.xlabel(join([xlabel " " xunits]), fontsize=labelsize)
	plt.ylabel(join([ylabel " " yunits]), fontsize=labelsize)
	ax = plt.gca()
	plt.setp(ax[:get_xticklabels](), fontsize=ticksize)
	plt.setp(ax[:get_yticklabels](), fontsize=ticksize)
	if (name == "NULL")
		plt.show()
	else
		plt.savefig(name,dpi=dpi)
		plt.close()
	end
	return im
end

function SeisPlot(in,extent::Extent;style="color",wiggle_fill_color="k",wiggle_line_color="k",wiggle_trace_increment=1,xcur=1.2,cmap="Greys",aspect="auto",pclip=98,vmin="NULL",vmax="NULL",yunits=" ",dpi=100,wbox=6,hbox=6,name="NULL",interpolation="Hanning",titlesize=20,labelsize=15,ticksize=10,fignum="NULL",plot_type="TX",fmax=100)
	im = SeisPlot(in,style=style,wiggle_fill_color=wiggle_fill_color,wiggle_line_color=wiggle_line_color,wiggle_trace_increment=wiggle_trace_increment,xcur=xcur,cmap=cmap,aspect=aspect,pclip=pclip,vmin=vmin,vmax=vmax,title=extent.title,xlabel=extent.label2,xunits=join(["(",extent.unit2,")"]),ylabel=extent.label1,yunits=join(["(",extent.unit1,")"]),ox=extent.o2,dx=extent.d2,oy=extent.o1,dy=extent.d1,dpi=dpi,wbox=wbox,hbox=hbox,name=name,interpolation=interpolation,titlesize=titlesize,labelsize=labelsize,ticksize=ticksize,fignum=fignum,plot_type=plot_type,fmax=fmax)
	return im
end

function SeisPlot(in::ASCIIString;style="color",wiggle_fill_color="k",wiggle_line_color="k",wiggle_trace_increment=1,xcur=1.2,cmap="Greys",aspect="auto",pclip=98,vmin="NULL",vmax="NULL",title=" ",xlabel=" ",xunits=" ",ylabel=" ",yunits=" ",ox=0,dx=1,oy=0,dy=1,dpi=100,wbox=6,hbox=6,name="NULL",interpolation="Hanning",titlesize=20,labelsize=15,ticksize=10,fignum="NULL",plot_type="TX",fmax=100)

	d,h = SeisRead(in)
	SeisPlot(d,style=style,wiggle_fill_color=wiggle_fill_color,wiggle_line_color=wiggle_line_color,wiggle_trace_increment=wiggle_trace_increment,xcur=xcur,cmap=cmap,aspect=aspect,pclip=pclip,vmin=vmin,vmax=vmax,title=title,xlabel=xlabel,xunits=xunits,ylabel=ylabel,yunits=yunits,ox=ox,dx=dx,oy=oy,dy=dy,dpi=dpi,wbox=wbox,hbox=hbox,name=name,interpolation=interpolation,titlesize=titlesize,labelsize=labelsize,ticksize=ticksize,fignum=fignum,plot_type=plot_type,fmax=fmax)

end
