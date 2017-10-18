# ReadWrite
Seismic.jl provides seismic data reading, writing and plotting. Currently, The conversions between SEIS data format and SEGY, SU and madagascar data format, bellow we give a simple example about converting SEGY data format to our internal data format.

## Example
```@example
using PyPlot, Seismic,Compat
download("http://seismic.physics.ualberta.ca/data/gom_cdp_nmo.su","gom_cdp_nmo.su");
SegyToSeis("gom_cdp_nmo.su","gom_cdp_nmo",format="su",input_type="ieee",swap_bytes=true);
d, h, ext = SeisRead("gom_cdp_nmo");
SeisPlot(d[:, 1:1000], cmap="PuOr", wbox=9)
savefig("usgs.svg"); nothing # hide
```
In the above example, we first download the data, then convert the data from SU data format to SEIS format, finally the data are plotted.
![](usgs.svg)
