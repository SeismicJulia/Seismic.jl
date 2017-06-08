
<a id='Plotting-1'></a>

# Plotting


Tools for figure plotting


<a id='SeisPlot-1'></a>

## SeisPlot

<a id='Seismic.SeisPlot' href='#Seismic.SeisPlot'>#</a>
**`Seismic.SeisPlot`** &mdash; *Function*.



```
SeisPlot(d [, extent]; <keyword arguments>)
```

Plot time-space, frequency-wavenumber or amplitude-frequency 2D seismic data `d` with color, wiggles or overlay.

**Arguments**

  * `d`: 2D data to plot.
  * `extent`: extent of the data (optional).

**Keyword arguments**

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

**Example**

```julia
julia> d, extent = SeisLinearEvents(); SeisPlot(d);
```

Credits: Aaron Stanton, 2015


<a target='_blank' href='https://github.com/fercarozzi/myseismicjulia/tree/c2832f8331d8b4cba573c54c2dd0183c518801d7/src/Plotting/SeisPlot.jl#L1-L55' class='documenter-source'>source</a><br>


<a id='SeisPlotCoordinates-1'></a>

## SeisPlotCoordinates

<a id='Seismic.SeisPlotCoordinates' href='#Seismic.SeisPlotCoordinates'>#</a>
**`Seismic.SeisPlotCoordinates`** &mdash; *Function*.



```
SeisPlotCoordinates(headers;<keyword arguments>)
```

Plot shot and receiver coordinates or fold map. Optionally can print to a file.

**Arguments**

  * `headers` :a vector where each element is of type Header

**Keyword arguments**

  * `style="sxsygxgy"` or "fold"
  * `cmap="Greys"` or "PuOr"
  * `aspect="auto"`: color image aspect ratio.
  * `pclip=98`: percentile for determining clip.
  * `vmin="NULL"`: minimum value used in colormapping data.
  * `vmax="NULL"`: maximum value used in colormapping data.
  * `title=" "`: title of plot.
  * `xlabel=" "`: label on x-axis.
  * `xunits=" "`: units of y-axis.
  * `ylabel=" "`: label on x-axis.
  * `yunits=" "`: units of y-axis.
  * `ox=0`: first point of x-axis.
  * `dx=1`: increment of x-axis.
  * `oy=0`: first point of y-axis.
  * `dy=1`: increment of y-axis.
  * `dpi=100`: dots-per-inch of figure.
  * `wbox=6`: width of figure in inches.
  * `hbox=6`: height of figure in inches.
  * `name="NULL"`: name of the figure to save (only if `name` is given).
  * `interpolation="none"`: interpolation method for colormapping data.
  * `titlesize=16`: size of title.
  * `labelsize=14`: size of labels on axis.
  * `ticksize=11`: size of labels on ticks.
  * `fignum="NULL"`: number of figure.

*Credits: Aaron Stanton, 2015*


<a target='_blank' href='https://github.com/fercarozzi/myseismicjulia/tree/c2832f8331d8b4cba573c54c2dd0183c518801d7/src/Plotting/SeisPlotCoordinates.jl#L1-L38' class='documenter-source'>source</a><br>

