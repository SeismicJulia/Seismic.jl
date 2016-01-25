using Seismic

download("http://certmapper.cr.usgs.gov/nersl/NPRA/seismic/1979/616_79/PROCESSED/616_79_PR.SGY","616_79_PR.SGY")
SegyToSeis("616_79_PR.SGY","616_79_PR.seis")
d,h,e = SeisRead("616_79_PR.seis")
imx = Seismic.ExtractHeader(h,"imx")
println("size(d)=",size(d))
println("imx=",imx)
SeisWrite("616_79_PR_b",d,h,e)

