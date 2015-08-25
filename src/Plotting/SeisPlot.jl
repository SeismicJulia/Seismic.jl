@doc """
Plot 2D data.
""" ->
function SeisPlot(in,param=Dict())

	canvas = get(param, "canvas", "NULL") # called by GUI, or REPL
	style = get(param,"style","color") # color or wiggles
	wiggle_fill_color = get(param,"wiggle_fill_color","k") # fill color for wiggles
	wiggle_line_color = get(param,"wiggle_line_color","k") # line color for wiggles
	wiggle_trace_increment = get(param,"wiggle_trace_increment",1) # can skip traces for wiggles
	xcur = get(param,"xcur",1.2) # wiggle excursion if wiggles selected
	cmap = get(param,"cmap","RdGy")
	aspect = get(param,"aspect","auto")
	pclip = get(param,"pclip",95)/100.
	vmin = get(param,"vmin",-pclip*maximum(abs(in[:])))
	vmax = get(param,"vmax",pclip*maximum(abs(in[:])))
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
	wbox = get(param,"wbox",4)
	hbox = get(param,"hbox",4)
	name = get(param,"name","NULL")
	interpolation = get(param,"interpolation","none")

	if (canvas == "NULL")
		plt.ion()
		if (!haskey(param,"fignum"))
			fig = plt.figure(figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
		else
			fig = plt.figure(num=get(param,"fignum",1), figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
		end
	else
		fig = canvas[:get_figure]()
	end

	# runs if we want a colour plot, imshow() draws desired plot with given parameters
	if (style != "wiggles")
		# when called by the GUI as opposed to the REPL
		if (canvas != "NULL")
			canvas[:imshow](in,cmap=cmap,vmin=vmin,vmax=vmax,extent=[ox,ox + (size(in,2)-1)*dx,oy + (size(in,1)-1)*dy,oy],aspect=aspect)
		else
			plt.imshow(in,cmap=cmap,vmin=vmin,vmax=vmax,extent=[ox,ox + (size(in,2)-1)*dx,oy + (size(in,1)-1)*dy,oy],aspect=aspect)
		end
	end

	# runs when we want a wiggle plot
	if (style != "color")
		# creates arrays x and y that are [0:xmax] and [y:ymax]
		y = oy+dy*[0:1:size(in,1)-1]
		x = ox+dx*[0:1:size(in,2)-1]
		# related to ???
		delta = wiggle_trace_increment*dx
		# set maximums
		hmin = minimum(x)
		hmax = maximum(x)
		dmax = maximum(in[:])
		# related to ???
		alpha = xcur*delta/dmax

		# loops x times and draws, and fills in, each trace
		for k = 1:wiggle_trace_increment:size(in,2)
			# arrays to help draw vertical lines individually
			x_vert = Float64[]
			y_vert = Float64[]
			# an array that contains the current trace's data
			s = in[:,k]
			# set the beginning and end to 0 to help the plot() function
			s[1]=0
			s[end]=0
			# an array with only the positive numbers from s, negative numbers returned to 0
			# helps with filling
			sp = (s+abs(s))/2

			# when called by the GUI
			if (canvas != "NULL")
				# draws the wiggle line of current trace, adds it to the plot
				wigline = lines.Line2D(s*alpha + x[k], y, color=wiggle_line_color)
				canvas[:add_line](wigline)
				# draws the vertical line through the current trace, adds it to the plot
				for i = 1:1:size(in,1)
					push!(x_vert, k)
					push!(y_vert, i)
				end
				vline = lines.Line2D(x_vert, y_vert, color=wiggle_line_color)
				canvas[:add_line](vline)
			else
				# all of the above work is done by plot() when called by pyplot in the REPL
				plt.plot( s*alpha + x[k],y,wiggle_line_color)
			end

			# unless overlay is called, the wiggle plot is filled
			if (style != "overlay") 
				if (canvas != "NULL")
					# fills between positive data and the vertical line of the trace
					canvas[:fill_between](sp*alpha + x[k],y_vert, color=wiggle_fill_color)
				else
					# same as above, for the REPL
					plt.fill(sp*alpha + x[k],y,wiggle_fill_color)
				end
			end
		end

		# set the visual parameters, axis markers, etc
		if (canvas != "NULL")
			canvas[:axis]([ox,ox + (size(in,2)-1)*dx,oy + (size(in,1)-1)*dy,oy])
		else
			plt.axis([ox,ox + (size(in,2)-1)*dx,oy + (size(in,1)-1)*dy,oy])
		end
	end

	if (canvas == "NULL")
		plt.title(title)
		plt.xlabel(join([xlabel " " xunits]))
		plt.ylabel(join([ylabel " " yunits]))
	end

	if (name == "NULL")
		plt.show()
	else  
		plt.savefig(name,dpi=dpi)
		plt.close()
	end
end
