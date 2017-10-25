<a name="logo"/>
<div align="center">
<a href="http://saig.physics.ualberta.ca/" target="_blank">
<img src="https://saig.physics.ualberta.ca/lib/tpl/dokuwiki/images/logo.png" alt="SAIG Logo" width="240" height="106"></img>
</a>
</div>

# Seismic.jl

[![Build Status](https://travis-ci.org/SeismicJulia/Seismic.jl.svg?branch=master)](https://travis-ci.org/SeismicJulia/Seismic.jl)

This module provides tools to read, write, and process 
seismic data. 
(The documentation can be found [here](http://seismicjulia.github.io/Seismic.jl)
and [here]
(http://seismic.physics.ualberta.ca).). At the moment, it is updated and tested against Julia 0.6

## Installation
To use this package you must first install the [Julia](http://julialang.org/downloads/) programming language. Once you have Julia you can install the Seismic package by typing ```Pkg.add("Seismic")``` on the Julia command line and then run ```Pkg.checkout("Seismic")``` to stay updated to the last version in this repository. 

## Basic usage
Once you have installed the package you can type `using Seismic` to start using
the functions. For example

```Julia
using PyPlot, Seismic;
download("http://seismic.physics.ualberta.ca/data/gom_cdp_nmo.su","gom_cdp_nmo.su");
SegyToSeis("gom_cdp_nmo.su","gom_cdp_nmo", format="su", input_type="ieee",swap_bytes=true);
d, h, e = SeisRead("gom_cdp_nmo");
SeisPlot(d[1:500, :]);
```
will produce this figure:

![plot1](http://seismic.physics.ualberta.ca/figures/616_79_PR.png)

## For developers: contributing to the package
If you want to fork the repository and contribute to Seismic.jl:
* New at GitHub? These [basic commands](http://seismic.physics.ualberta.ca/docs/git_basic_commands.pdf) 
and this [dictionary](http://seismic.physics.ualberta.ca/docs/git_dictionary.pdf) might help.
* This [tutorial](http://seismic.physics.ualberta.ca/docs/develop_SeismicJulia.pdf) provides the basics 
steps you need to follow in order to fork the main repository, change the source code in your forked 
repository, commit the changes, and make pull requests using GitHub.
* For contributions to the package, please follow the general guidelines given here: 
[Modifications.md](https://github.com/SeismicJulia/Seismic.jl/blob/master/Modifications.md).
