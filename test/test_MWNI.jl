using Compat
@compat using Seismic
using Base.Test

d,e = SeisLinearEvents()
ddec = SeisDecimate(d,perc=50)

dmwni = SeisMWNI(ddec,dt=0.004,fmax=100,Niter_external=5,Niter_internal=10)

# test that quality factor is greater than 10 Decibels
quality_factor = 10*log10(norm(d[:],2)/norm(dmwni[:]-d[:],2))
println("Quality factor = ",quality_factor)
@test quality_factor > 5
