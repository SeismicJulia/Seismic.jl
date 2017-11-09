# Seismic.jl

Seismic data processing, imaging and plotting

---

## Overview

Seismic.jl provides tools to **process**, **image**, and **plot** reflection seismic data in the Julia language.

### _Convert data to a simple format_

Data and headers are stored separately as _filename@data@_ and _filename@headers@_ for 
simplicity. Functions are available to convert from popular formats such as SEGY, SU and RSF.

### _Manipulate data_

Contains functions for geometry calculation, sorting, windowing, patching/un-patching, and processing keyed on header word.

### _Plot publication quality figures_

Produce color and wiggle plots using PyPlot.jl

---

## Installation

To install Seismic.jl type:

```no-highlight
Pkg.add("Seismic")
```
```no-highlight
Pkg.checkout("Seismic")
```
---

## Getting started

To start using the functions simply type `using Seismic`. Below is a simple demonstration of
the plotting functionality: 

```no-highlight
using PyPlot,Seismic

param = Dict(:nt=>500, :nx1=>500, :tau=>[0.4, 1.0], :p1=>[-.00003, 0.00008], :amp=>[1., -0.5], :f0=>20.0)

d,ext = SeisLinearEvents(;param...)

plotpar = Dict(:style=>"overlay",
           :wiggle_trace_increment=>10,
           :xcur=>0.8,
           :aspect=>"auto",
           :xlabel=>"X",:xunits=>"meters",:ox=>0,:dx=>10,
           :ylabel=>"Time",:yunits=>"seconds",:oy=>0,:dy=>0.004,
           :wbox=>8,:hbox=>5,
           :cmap=>"seismic");

plotpar[:style]="color"; plotpar[:title]="color"; plotpar[:name]="plot1"; 
SeisPlot(d;plotpar...);
plotpar[:style]="wiggles"; plotpar[:title]="wiggles"; plotpar[:name]="plot2"; 
SeisPlot(d;plotpar...);
plotpar[:style]="overlay"; plotpar[:title]="overlay"; plotpar[:name]="plot3"; 
SeisPlot(d;plotpar...);
```
will produce these three `.eps` files:

![plot1](http://www.ualberta.ca/~kstanton/files/plot1.png "color")

![plot2](http://www.ualberta.ca/~kstanton/files/plot2.png "wiggles")

![plot3](http://www.ualberta.ca/~kstanton/files/plot3.png "overlay")

