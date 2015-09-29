function SeisMap(headers,param)

	style = get(param,"style","sxsygxgy") # type of map
	fignum = get(param,"fignum",1)
	cmap = get(param,"cmap","RdGy") # color map for fold plot
	aspect = get(param,"aspect","auto")
	vmin = get(param,"vmin",-1) # min clip value for fold plot
	vmax = get(param,"vmax",1) # max clip value for fold plot
	title  = get(param,"title"," ")
	xlabel = get(param,"xlabel"," ")
	xunits = get(param,"xunits"," ")
	ylabel = get(param,"ylabel"," ")
	yunits = get(param,"yunits"," ")
	ox = get(param,"ox",0)
	dx = get(param,"dx",1)
	oy = get(param,"oy",0)
	dy = get(param,"dy",1)
	dpi = get(param,"dpi",100)
	wbox = get(param,"wbox",8)
	hbox = get(param,"hbox",8)
	name = get(param,"name","NULL")
	interpolation = get(param,"interpolation","none")

	fig = plt.figure(num=fignum, figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
	if (style == "sxsygxgy")
		ntrace = length(headers)
		sx = zeros(ntrace)
		sy = zeros(ntrace)
		gx = zeros(ntrace)
		gy = zeros(ntrace)
		for itrace = 1 : ntrace
			sx[itrace] = headers[itrace].sx;
			sy[itrace] = headers[itrace].sy;
			gx[itrace] = headers[itrace].gx;
			gy[itrace] = headers[itrace].gy;
		end
		plot(gx,gy,linestyle="None",marker="^",markersize=5,color="b");
		plot(sx,sy,linestyle="None",marker="*",markersize=8,color="r");
		#plt.axis([ox,ox + (size(in,2)-1)*dx,oy + (size(in,1)-1)*dy,oy])
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

end
