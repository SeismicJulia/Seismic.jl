# Seismic.jl

Seismic data processing, imaging and plotting

---

## Overview

Seismic.jl provides tools to **process**, **image**, and **plot** reflection seismic data in the Julia language.

### Read and write SEGY, SU, and RSF files

Convert to and from SEGY, SU and RSF file formats

### Simple file format

Data and headers are stored separately as _filename.seisd_ and _filename.seish_ for simplicity.

### Data manipulation

Functions for geometry calculation, sorting, windowing, patching/un-patching, and processing keyed on header word.

### Plotting 

Color and wiggle plots using PyPlot.jl

---

## Installation

To install Seismic.jl type:

```no-highlight
Pkg.add("Seismic")
```

---

## Getting started

To start using the functions simply type `using Seismic`. Below is a simple demonstration of
the plotting functionality: 

```no-highlight
using PyPlot,Seismic

param = ["nt"=>500,"nx1"=>500,
	 "tau1"=>[0.4 1.0],"tau2"=>[0. 0.],"tau3"=>[0. 0.],"tau4"=>[0. 0.],
	 "v1"=>[3500. -4000],"v2"=>[99999. 99999.],"v3"=>[99999. 99999.],"v4"=>[99999. 99999.],
     "amp"=>[1. -0.5], "f0"=>[20. 20.]];
d = SeisLinearEvents(param);

plotpar = {"style"=>"overlay",
           "wiggle_trace_increment"=>10,
           "xcur"=>0.8,
           "vmin"=>-2,"vmax"=>2,
           "aspect"=>"auto",
           "xlabel"=>"X","xunits"=>"meters","ox"=>0,"dx"=>10,
           "ylabel"=>"Time","yunits"=>"seconds","oy"=>0,"dy"=>0.004,
           "wbox"=>8,"hbox"=>5,
           "cmap"=>"seismic"};

plotpar["style"]="color"; plotpar["title"]="color"; plotpar["name"]="plot1"; 
SeisPlot(d[:,:],plotpar);
plotpar["style"]="wiggles"; plotpar["title"]="wiggles"; plotpar["name"]="plot2"; 
SeisPlot(d[:,:],plotpar);
plotpar["style"]="overlay"; plotpar["title"]="overlay"; plotpar["name"]="plot3"; 
SeisPlot(d[:,:],plotpar);
```
will produce these three `.eps` files:

![plot1](http://www.ualberta.ca/~kstanton/files/plot1.png "color")

![plot2](http://www.ualberta.ca/~kstanton/files/plot2.png "wiggles")

![plot3](http://www.ualberta.ca/~kstanton/files/plot3.png "overlay")

