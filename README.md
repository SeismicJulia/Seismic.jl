<a name="logo"/>
<div align="center">
<a href="http://saig.physics.ualberta.ca/" target="_blank">
<img src="https://saig.physics.ualberta.ca/lib/tpl/dokuwiki/images/logo.png" alt="SAIG Logo" width="240" height="106"></img>
</a>
</div>

# Seismic.jl

[![Build Status](https://travis-ci.org/SeismicJulia/Seismic.jl.svg?branch=master)](https://travis-ci.org/SeismicJulia/Seismic.jl)

This module provides tools to read, write, process, and plot 3D reflection 
seismic data. 
[//]: # (The documentation can be found [here](http://seismic.physics.ualberta.ca).)

## Installation
To use this package you must first install [the Julia programming language](http://julialang.org). Once you have Julia you can download and install the Seismic package by typing ```Pkg.add("Seismic")``` on the Julia command line.

## Basic usage
Once you have installed the package you can type `using Seismic` to start using
the functions. For example

```Julia
using PyPlot, Seismic;
download("http://certmapper.cr.usgs.gov/nersl/NPRA/seismic/1979/616_79/PROCESSED/616_79_PR.SGY", "616_79_PR.SGY");
SegyToSeis("616_79_PR.SGY", "616_79_PR.seis");
d, h, e = SeisRead("616_79_PR.seis");
SeisPlot(d[1:500, :], e, cmap="PuOr", wbox=9);
```
will produce this figure:

![plot1](http://seismic.physics.ualberta.ca/figures/616_79_PR.png)

## For developers: contributing to the package
If you want to fork the repository and contribute to Seismic.jl: 
* Follow [This tutorial](http://seismic.physics.ualberta.ca/docs/develop_SeismicJulia.pdf) that explains the basics of how to do it.; and
* Follow the general guidelines described here: [Modifications.md](https://github.com/SeismicJulia/Seismic.jl/blob/master/Modifications.md).
