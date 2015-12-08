<a name="logo"/>
<div align="center">
<a href="http://saig.physics.ualberta.ca/" target="_blank">
<img src="https://saig.physics.ualberta.ca/lib/tpl/dokuwiki/images/logo.png" alt="SAIG Logo" width="240" height="106"></img>
</a>
</div>

# Seismic.jl

[![Build Status](https://travis-ci.org/SeismicJulia/Seismic.jl.svg?branch=master)](https://travis-ci.org/SeismicJulia/Seismic.jl)

This module provides tools to read, write, process, and plot 3D reflection 
seismic data. The documentation can be found [here](http://seismic.physics.ualberta.ca).

## Installation
To use this package you must first install [the Julia programming language](http://julialang.org). Once you have Julia you can download and install the Seismic package by typing ```Pkg.add("Seismic")``` on the Julia command line.

## Basic usage
Once you have installed the package you can type `using Seismic` to start using
the functions. For example

```
using PyPlot,Seismic

param = ["nt"=>500,"nx1"=>500,
	 "tau1"=>[0.4 1.0],"tau2"=>[0. 0.],"tau3"=>[0. 0.],"tau4"=>[0. 0.],
	 "v1"=>[3500. -4000],"v2"=>[99999. 99999.],"v3"=>[99999. 99999.],"v4"=>[99999. 99999.],
     "amp"=>[1. -0.5], "f0"=>[20. 20.]];
d = SeisLinearEvents(param);

plotpar = ["style"=>"overlay",
           "wiggle_trace_increment"=>10,
           "xcur"=>0.8,
           "vmin"=>-2,"vmax"=>2,
           "aspect"=>"auto",
           "xlabel"=>"X","xunits"=>"meters","ox"=>0,"dx"=>10,
           "ylabel"=>"Time","yunits"=>"seconds","oy"=>0,"dy"=>0.004,
           "wbox"=>8,"hbox"=>5,
           "cmap"=>"seismic"];

plotpar["style"]="color"; plotpar["title"]="color"; plotpar["name"]="plot1"; 
SeisPlot(d[:,:],plotpar);
plotpar["style"]="wiggles"; plotpar["title"]="wiggles"; plotpar["name"]="plot2"; 
SeisPlot(d[:,:],plotpar);
plotpar["style"]="overlay"; plotpar["title"]="overlay"; plotpar["name"]="plot3"; 
SeisPlot(d[:,:],plotpar);
```
will produce these three `.png` files:

![plot1](http://www.ualberta.ca/~kstanton/files/plot1.png "color")

![plot2](http://www.ualberta.ca/~kstanton/files/plot2.png "wiggles")

![plot3](http://www.ualberta.ca/~kstanton/files/plot3.png "overlay")
