# Seismic.jl

Seismic signal analysis and imaging

---

## Overview

Seismic.jl provides tools to **process**, **image**, and **plot** reflection seismic data in the Julia language.

### _Convert data to a simple format_

Data and headers are stored separately as _filename.seisd_ and _filename.seish_ for
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

## Index
```@index
```
