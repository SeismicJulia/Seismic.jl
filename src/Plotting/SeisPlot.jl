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
* interpolation="none"
* titlesize=20
* labelsize=15
* ticksize=10
* fignum="NULL"

**OUT**  

*Credits: Aaron Stanton, 2015*

"""
function SeisPlot(in;style="color",wiggle_fill_color="k",wiggle_line_color="k",wiggle_trace_increment=1,xcur=1.2,cmap="Greys",aspect="auto",pclip=98,vmin="NULL",vmax="NULL",title=" ",xlabel=" ",xunits=" ",ylabel=" ",yunits=" ",ox=0,dx=1,oy=0,dy=1,dpi=100,wbox=6,hbox=6,name="NULL",interpolation="none",titlesize=20,labelsize=15,ticksize=10,fignum="NULL")

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
	if (style != "wiggles")
		im = plt.imshow(in,cmap=cmap,vmin=a,vmax=b,extent=[ox,ox + (size(in,2)-1)*dx,oy + (size(in,1)-1)*dy,oy],aspect=aspect)
	end
	if (style != "color")
		y = oy+dy*[0:1:size(in,1)-1]
		x = ox+dx*[0:1:size(in,2)-1]
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

function SeisPlot(in::ASCIIString;style="color",wiggle_fill_color="k",wiggle_line_color="k",wiggle_trace_increment=1,xcur=1.2,cmap="Greys",aspect="auto",pclip=98,vmin="NULL",vmax="NULL",title=" ",xlabel=" ",xunits=" ",ylabel=" ",yunits=" ",ox=0,dx=1,oy=0,dy=1,dpi=100,wbox=6,hbox=6,name="NULL",interpolation="none",titlesize=20,labelsize=15,ticksize=10,fignum="NULL")

	d,h = SeisRead(in)
	SeisPlot(d,style=style,wiggle_fill_color=wiggle_fill_color,wiggle_line_color=wiggle_line_color,wiggle_trace_increment=wiggle_trace_increment,xcur=xcur,cmap=cmap,aspect=aspect,pclip=pclip,vmin=vmin,vmax=vmax,title=title,xlabel=xlabel,xunits=xunits,ylabel=ylabel,yunits=yunits,ox=ox,dx=dx,oy=oy,dy=dy,dpi=dpi,wbox=wbox,hbox=hbox,name=name,interpolation=interpolation,titlesize=titlesize,labelsize=labelsize,ticksize=ticksize,fignum=fignum)

end
