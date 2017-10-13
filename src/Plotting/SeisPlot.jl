"""
    SeisPlot(d [, extent]; <keyword arguments>)

Plot time-space, frequency-wavenumber or amplitude-frequency 2D seismic data `d`
with color, wiggles or overlay.

# Arguments
* `d`: 2D data to plot.
* `extent`: extent of the data (optional).

# Keyword arguments
* `plot_type="TX"`: `"TX"` for time-space plot; `"FK"` for frequency-wavenumber
plot or `"Amplitude"` for amplitude spectrum plot.
* `style="color"`: style of the plot: `"color"`, `"wiggles"` or `"overlay"`.
* `cmap="PuOr"`: colormap for  `"color"` or `"overlay"` style.
* `pclip=98`: percentile for determining clip.
* `vmin="NULL"`: minimum value used in colormapping data.
* `vmax="NULL"`: maximum value used in colormapping data.
* `aspect="auto"`: color image aspect ratio.
* `interpolation="Hanning"`: interpolation method for colormapping data.
* `fmax=100`: maximum frequency for `"FK"` or `"Amplitude"` plot.
* `wiggle_fill_color="k"`: color for filling the positive wiggles.
* `wiggle_line_color="k"`: color for wiggles' lines.
* `wiggle_trace_increment=1`: increment for wiggle traces.
* `xcur=1.2`: wiggle excursion in traces corresponding to clip.
* `scal="NULL"`: scale for wiggles.
* `title=" "`: title of plot.
* `titlesize=16`: size of title.
* `xlabel=" "`: label on x-axis.
* `xunits=" "`: units of y-axis.
* `ylabel=" "`: label on x-axis.
* `yunits=" "`: units of y-axis.
* `labelsize=14`: size of labels on axis.
* `ox=0`: first point of x-axis.
* `dx=1`: increment of x-axis.
* `oy=0`: first point of y-axis.
* `dy=1`: increment of y-axis.
* `xticks="NULL"`: ticks on x-axis.
* `yticks="NULL"`: ticks on y-axis.
* `xticklabels="NULL"`: labels on ticks of x-axis.
* `yticklabels="NULL"`: labels on ticks of y-axis.
* `ticksize=11`: size of labels on ticks.
* `fignum="NULL"`: number of figure.
* `wbox=6`: width of figure in inches.
* `hbox=6`: height of figure in inches.
* `dpi=100`: dots-per-inch of figure.
* `name="NULL"`: name of the figure to save (only if `name` is given).

# Example
```julia
julia> d, extent = SeisLinearEvents(); SeisPlot(d);
```

Credits: Aaron Stanton, 2015
"""
function SeisPlot{T<:Real}(d::Array{T,2}; plot_type="TX", style="color",
                           cmap="PuOr", pclip=98, vmin="NULL", vmax="NULL",
                           aspect="auto", interpolation="Hanning", fmax=100,
                           wiggle_fill_color="k", wiggle_line_color="k",
                           wiggle_trace_increment=1, xcur=1.2, scal="NULL",
                           title=" ", titlesize=16, xlabel=" ", xunits=" ",
                           ylabel=" ", yunits=" ", labelsize=14, ox=0, dx=1,
                           oy=0, dy=1, xticks="NULL", yticks="NULL",
                           xticklabels="NULL", yticklabels="NULL", ticksize=11,
                           fignum="NULL", wbox=6, hbox=6, dpi=100, name="NULL")
    if (vmin=="NULL" || vmax=="NULL")
        if (pclip<=100)
	    a = -quantile(abs.(d[:]), (pclip/100))
	else
	    a = -quantile(abs.(d[:]), 1)*pclip/100
	end
	b = -a
    else
	a = vmin
	b = vmax
    end
    plt.ion()
    if (fignum == "NULL")
	fig = plt.figure(figsize=(wbox, hbox), dpi=dpi, facecolor="w",
                           edgecolor="k")
    else
	fig = plt.figure(num=fignum, figsize=(wbox, hbox), dpi=dpi,
                           facecolor="w", edgecolor="k")
    end
    if plot_type == "TX"
	if (style != "wiggles")
	    im = plt.imshow(d, cmap=cmap, vmin=a, vmax=b,
                            extent=[ox - dx/2,ox + (size(d,2)-1)*dx + dx/2,
                                    oy + (size(d,1)-1)*dy,oy],
                            aspect=aspect, interpolation=interpolation)
	end
	if (style != "color")
            style=="wiggles" ? margin = dx : margin = dx/2
	    y = oy+dy*collect(0:1:size(d, 1)-1)
	    x = ox+dx*collect(0:1:size(d, 2)-1)
	    delta = wiggle_trace_increment*dx
	    hmin = minimum(x)
	    hmax = maximum(x)
            dmax = maximum(abs.(d[:]))
	    alpha = xcur*delta
            scal=="NULL" ? alpha = alpha/maximum(abs.(d[:])) : alpha=alpha*scal
	    for k = 1:wiggle_trace_increment:size(d, 2)
		x_vert = Float64[]
		y_vert = Float64[]
		sc = x[k] * ones(size(d, 1))
    s  = d[:,k]*alpha + sc
		im = plt.plot( s, y, wiggle_line_color)
		if (style != "overlay")
		    plt.fill_betweenx(y, sc, s, where=s.>sc, facecolor=wiggle_line_color)
		end
	    end
	    plt.axis([ox - margin, ox + (size(d, 2)-1)*dx + margin,
                      oy + (size(d, 1)-1)*dy, oy])
	end
    elseif plot_type == "FK"
	xlabel = "Wavenumber"
	xunits = "(1/m)"
	ylabel = "Frequency"
	yunits = "Hz"
	dk = 1/dx/size(d[:,:], 2)
	kmin = -dk*size(d[:,:], 2)/2
	kmax =  dk*size(d[:,:], 2)/2
	df = 1/dy/size(d[:,:], 1)
	FMAX = df*size(d[:,:], 1)/2
	if fmax > FMAX
	    fmax = FMAX
	end
	nf = convert(Int32, floor((size(d[:, :], 1)/2)*fmax/FMAX))
	D = abs.(fftshift(fft(d[:, :])))
	D = D[round(Int,end/2):round(Int,end/2)+nf, :]
	if (vmin=="NULL" || vmax=="NULL")
	    a = 0.
	    if (pclip<=100)
		b = quantile(abs.(D[:]), (pclip/100))
	    else
		b = quantile(abs.(D[:]), 1)*pclip/100
	    end
	end
	im = plt.imshow(D, cmap=cmap, vmin=a, vmax=b, extent=[kmin,kmax,fmax,0],
                        aspect=aspect, interpolation=interpolation)
    elseif plot_type == "Amplitude"
	xlabel = "Frequency"
	xunits = "(Hz)"
	ylabel = "Amplitude"
	yunits = ""
	nx = size(d[:,:], 2)
	df = 1/dy/size(d[:, :], 1)
	FMAX = df*size(d[:, :], 1)/2
	if fmax > FMAX
	    fmax = FMAX
	end
	nf = convert(Int32, floor((size(d[:, :], 1)/2)*fmax/FMAX))
	y = fftshift(sum(abs.(fft(d[:, :], 1)), 2))/nx
	y = y[round(Int,end/2):round(Int, end/2)+nf]
	norm = maximum(y[:])
	if (norm > 0.)
		y = y/norm
	end
	x = collect(0:df:fmax)
	im = plt.plot(x, y)
	plt.title(title)
	plt.xlabel(join([xlabel " " xunits]))
	plt.ylabel(join([ylabel " " yunits]))
	plt.axis([0, fmax, 0, 1.1])
    else
	error("plot_type not recognized.")
    end
    plt.title(title, fontsize=titlesize)
    plt.xlabel(join([xlabel " " xunits]), fontsize=labelsize)
    plt.ylabel(join([ylabel " " yunits]), fontsize=labelsize)
    xticks == "NULL" ? nothing : plt.xticks(xticks)
    yticks == "NULL" ? nothing : plt.yticks(yticks)
    ax = plt.gca()
    xticklabels == "NULL" ? nothing : ax[:set_xticklabels](xticklabels)
    yticklabels == "NULL" ? nothing : ax[:set_yticklabels](yticklabels)
    plt.setp(ax[:get_xticklabels](), fontsize=ticksize)
    plt.setp(ax[:get_yticklabels](), fontsize=ticksize)
    if (name == "NULL")
	plt.show()
    else
	plt.savefig(name, dpi=dpi)
	plt.close()
    end
    return im
end

function SeisPlot{T<:Real}(d::Array{T,2}, extent::Extent;
                           plot_type="TX", style="color",
                           cmap="PuOr", pclip=98, vmin="NULL", vmax="NULL",
                           aspect="auto", interpolation="Hanning", fmax=100,
                           wiggle_fill_color="k", wiggle_line_color="k",
                           wiggle_trace_increment=1, xcur=1.2, scal="NULL",
                           title=" ", titlesize=16, xlabel=" ", xunits=" ",
                           ylabel=" ", yunits=" ", labelsize=14, ox=0, dx=1,
                           oy=0, dy=1, xticks="NULL", yticks="NULL",
                           xticklabels="NULL", yticklabels="NULL", ticksize=11,
                           fignum="NULL", wbox=6, hbox=6, dpi=100, name="NULL")
    im = SeisPlot(d, plot_type=plot_type, style=style, cmap=cmap, pclip=pclip,
                  vmin=vmin, vmax=vmax, aspect=aspect,
                  interpolation=interpolation, fmax=fmax,
                  wiggle_fill_color=wiggle_fill_color,
                  wiggle_line_color=wiggle_line_color,
                  wiggle_trace_increment=wiggle_trace_increment, xcur=xcur,
                  scal=scal, title=extent.title,  titlesize=titlesize,
                  xlabel=extent.label2, xunits=join(["(",extent.unit2,")"]),
                  ylabel=extent.label1, yunits=join(["(",extent.unit1,")"]),
                  labelsize=labelsize, ox=extent.o2, dx=extent.d2, oy=extent.o1,
                  dy=extent.d1, xticks=xticks, yticks=yticks,
                  xticklabels=xticklabels, yticklabels=yticklabels,
                  ticksize=ticksize, wbox=wbox, hbox=hbox, dpi=dpi, name=name,
                  fignum=fignum)
    return im
end
