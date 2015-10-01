# Seismic.jl

Tools to process and image seismic data in the Julia language.

---

## Overview

Seismic.jl is a **fast**, **simple** and **downright gorgeous** static site
generator that's geared towards building project documentation. Documentation
source files are written in Markdown, and configured with a single YAML
configuration file.

### Read and write Segy, SU, and RSF files

Builds completely static HTML sites that you can host on GitHub pages, Amazon
S3, or anywhere else you choose.


---

## Installation

To install Seismic.jl type:

```bash
Pkg.add("Seismic")
```

---

## Getting started

To start using the functions simply type `using Seismic`. Below is a simple demonstration of
the plotting functionality: 

```bash
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

