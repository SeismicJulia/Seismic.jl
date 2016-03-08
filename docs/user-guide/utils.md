# Manipulating Seismic data disk files

Often seismic data processing requires us to manipulate large disk files. These manipulations
include windowing, sorting, and processing small portions of a large dataset sequentially or 
in parallel. Below are some functions to manipulate .seisd/.seish files

---

## SeisWindow

Windowing a .seisd/.seish file based on one or more keywords and/or the first dimension (time or depth).

## SeisSort

## SeisProcess

This function give the ability to apply a process to a .seisd/.seish file all at once, or in 
groups of n traces, or based on a change in one or more header keys. For example, applying NMO
correction to individual CMP gathers is achieved by:

```no-highlight
SeisProcess("foo","bar",{"group"=>"key","key"=>["imx"],"f"=>[SeisNMO],"tnmo"=>[0. 9999.],"vnmo"=>[1500 1500]]})
```
Here the parameter "f"=>[SeisNMO] refers to only one process to be applied, but more than one
process can be specified. 

## SeisPatch

This function breaks the data into a number of multidimensional patches for processing on overlapping 
windows in up to 5 dimensions.

## SeisUnPatch

This function reassembles patches of overlapping data windows into a single file. Tapering is applied to 
overlapping regions to ensure partition of unity.

## SeisPatchProcess

This functions wraps both the `SeisPatch`, `SeisProccess` and `SeisUnPatch` functions so parallel processing 
can be applied to overlapping multidimensional windows in a single call.

## SeisGeometry

This function acts on the .seish file to update geometry headers, rotates survey coordinates to a North-South East-West 
coordinate frame and bins header values into integer values. 

## SeisBin

This function creates a five dimensional regularly sampled volume of zeros and maps the input data file into this 
volume by nearest neighbour binning.

