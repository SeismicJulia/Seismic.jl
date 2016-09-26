"""
**SeisPlotCoordinates**

*Plot shot and receiver coordinates or fold map. Optionally can print to a file. *

**IN**   
* h: headers (a vector where each element is of type Header) 

* style="sxsygxgy" or "fold"
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
function SeisPlotCoordinates(headers;style="sxsygxgy",cmap="Greys",aspect="auto",pclip=98,vmin="NULL",vmax="NULL",title=" ",xlabel=" ",xunits=" ",ylabel=" ",yunits=" ",ox=0,dx=1,oy=0,dy=1,dpi=100,wbox=6,hbox=6,name="NULL",interpolation="none",titlesize=20,labelsize=15,ticksize=10,fignum="NULL")

	plt.ion()
	if (fignum == "NULL")
		fig = plt.figure(figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
	else
		fig = plt.figure(num=fignum,figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
	end
	if (style == "sxsygxgy")
		sx = Seismic.ExtractHeader(headers,"sx")
		sy = Seismic.ExtractHeader(headers,"sy")
		gx = Seismic.ExtractHeader(headers,"gx")
		gy = Seismic.ExtractHeader(headers,"gy")
		im = plt.scatter(gx,gy,marker="^",s=5,color="b");
		im = plt.scatter(sx,sy,marker="*",s=8,color="r");
		#plt.axis([ox,ox + (size(in,2)-1)*dx,oy + (size(in,1)-1)*dy,oy])
	elseif (style == "fold")
		imx = Seismic.ExtractHeader(headers,"imx")
		imy = Seismic.ExtractHeader(headers,"imy")
		imx_min = minimum(imx)
		imx_max = maximum(imx)
		nmx = imx_max - imx_min + 1
		imy_min = minimum(imy)
		imy_max = maximum(imy)
		nmy = imy_max - imy_min + 1
		fold = zeros(nmy,nmx)
		for itrace = 1 : length(imx)
			imx_index = imx[itrace] - imx_min + 1
			imy_index = imy[itrace] - imy_min + 1
			fold[imy_index,imx_index] += 1
		end
		if (vmin=="NULL" || vmax=="NULL")
			if (pclip<=100)
				a = -quantile(abs(fold[:]),(pclip/100))
			else
				a = -quantile(abs(fold[:]),1)*pclip/100
			end
			b = -a
		else
			a = vmin
			b = vmax
		end
		im = plt.imshow(fold,cmap=cmap,vmin=a,vmax=b,extent=[ox,ox + (size(fold,2)-1)*dx,oy + (size(fold,1)-1)*dy,oy],aspect=aspect)
		ax = plt.gca()
		ax[:invert_yaxis]();
		plt.colorbar()
	else
		println("no other plotting style defined yet")
	end
	
	plt.title(title)
	plt.xlabel(join([xlabel " " xunits]))
	plt.ylabel(join([ylabel " " yunits]))

	if (name == "NULL")
		plt.show()
	else  
		plt.savefig(name,dpi=dpi)
		plt.close()
	end
	return im
end
