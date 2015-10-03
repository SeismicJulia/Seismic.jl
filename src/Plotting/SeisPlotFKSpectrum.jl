function SeisPlotFKSpectrum(in,param=Dict())

	fignum = get(param,"fignum",1)
	cmap = get(param,"cmap","jet")
	aspect = get(param,"aspect","auto")
	fmax = get(param,"fmax",100) # maximum frequency to plot
	pclip = get(param,"pclip",95)/100.
	vmin = get(param,"vmin",0)
	vmax = get(param,"vmax",pclip*maximum(abs(fft(in[:,:]))))
	title  = get(param,"title"," ")
	xlabel = get(param,"xlabel","Wavenumber")
	xunits = get(param,"xunits","(1/m)")
	ylabel = get(param,"ylabel","Frequency")
	yunits = get(param,"yunits","(Hz)")
	ox = get(param,"ox",0)
	dx = get(param,"dx",10)
	oy = get(param,"oy",0)
	dy = get(param,"dy",0.001)
	dpi = get(param,"dpi",100)
	wbox = get(param,"wbox",6)
	hbox = get(param,"hbox",6)
	name = get(param,"name","NULL")
	interpolation = get(param,"interpolation","none")
	dk = 1/dx/size(in[:,:],2)
	kmin = -dk*size(in[:,:],2)/2
	kmax =  dk*size(in[:,:],2)/2
	df = 1/dy/size(in[:,:],1)
	FMAX = df*size(in[:,:],1)/2 
	nf = convert(Int32,floor((size(in[:,:],1)/2)*fmax/FMAX))
	D = abs(fftshift(fft(in[:,:])))
	D = D[int(end/2):int(end/2)+nf,:]

	plt.ion()
	if (!haskey(param,"fignum"))
		fig = plt.figure(figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
	else
		fig = plt.figure(num=get(param,"fignum",1), figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
	end

	plt.imshow(D,cmap=cmap,vmin=vmin,vmax=vmax,extent=[kmin,kmax,fmax,0],aspect=aspect)
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
