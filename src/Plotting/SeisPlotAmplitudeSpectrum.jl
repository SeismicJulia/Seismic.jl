function SeisPlotAmplitudeSpectrum(in,param=Dict())

	fignum = get(param,"fignum",1)
	aspect = get(param,"aspect",2)
	fmax = get(param,"fmax",100) # maximum frequency to plot
	title  = get(param,"title"," ")
	xlabel = get(param,"xlabel","Frequency")
	xunits = get(param,"xunits","(Hz)")
	ylabel = get(param,"ylabel","Amplitude")
	yunits = get(param,"yunits"," ")
	ot = get(param,"ot",0)
	dt = get(param,"dt",0.001)
	dpi = get(param,"dpi",100)
	wbox = get(param,"wbox",6)
	hbox = get(param,"hbox",6)
	name = get(param,"name","NULL")
	interpolation = get(param,"interpolation","none")
	nx = size(in[:,:],2)
	df = 1/dt/size(in[:,:],1)
	FMAX = df*size(in[:,:],1)/2 
	nf = convert(Int32,floor((size(in[:,:],1)/2)*fmax/FMAX))
	y = fftshift(sum(abs(fft(in[:,:],1)),2))/nx
	y = y[int(end/2):int(end/2)+nf]
	norm = maximum(y[:])
	if (norm > 0.)
		y = y/norm
	end
	x = [0:df:fmax]

	plt.ion()
	if (!haskey(param,"fignum"))
		fig = plt.figure(figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
	else
		fig = plt.figure(num=get(param,"fignum",1), figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
	end

	plt.plot(x,y)
	plt.title(title)
	plt.xlabel(join([xlabel " " xunits]))
	plt.ylabel(join([ylabel " " yunits]))
	plt.axis([0,fmax,0,1.1])

	if (name == "NULL")
		plt.show()
	else  
		plt.savefig(name,dpi=dpi)
		plt.close()
	end

end
