# ReadWrite
Seismic.jl provides seismic data reading, writing and plotting. Currently, The conversions between .seis data format and .segy, .su and madagascar data format, bellow we give a simple example about converting SEGY data format to our internal data format.

## Example
```@example
using PyPlot, Seismic,Compat
download("http://certmapper.cr.usgs.gov/nersl/NPRA/seismic/1979/616_79/PROCESSED/616_79_PR.SGY", "616_79_PR.SGY");
SegyToSeis("616_79_PR.SGY", "616_79_PR.seis");
d, h, e = SeisRead("616_79_PR.seis");
SeisPlot(d[1:500, :], e, cmap="PuOr", wbox=9);
savefig("usgs.svg"); nothing # hide
```
In above example, we first download the data from USGS's website, then convert the data from SEGY data format to our internal format, finally the data are plotted.
![](usgs.svg)
