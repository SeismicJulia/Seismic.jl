using Seismic
using Base.Test

param = {"nt"=>500,"dt"=>0.004,"nx1"=>500,
"tau1"=>[0.4 1.0],"tau2"=>[0. 0.],"tau3"=>[0. 0.],"tau4"=>[0. 0.],
"v1"=>[3500. -4000],"v2"=>[99999. 99999.],"v3"=>[99999. 99999.],"v4"=>[99999. 99999.],
"amp"=>[1. -0.5], "f0"=>[20. 20.]};
d = SeisLinearEvents(param);

nt = size(d,1); 
nx = size(d[:,:],2)
h = Header[]
for ix = 1:nx
	push!(h,Seismic.InitSeisHeader())
	h[ix].tracenum = ix;
	h[ix].d1 = 0.004;
	h[ix].n1 = nt;
	h[ix].imx = ix;
end
param = {"mode"=>"random","perc"=> 50}
ddec,h = SeisDecimate(d,h,param)

param = {"style"=>"mxmyhxhy","Niter"=>100,"fmax"=>60,"padt"=>2,"padx"=>2}
dpocs,h = SeisPOCS(ddec,h,param)

# test that quality factor is greater than 10 Decibels
quality_factor = 10*log10(norm(d[:],2)/norm(dpocs[:]-d[:],2))
println("Quality factor = ",quality_factor)
@test quality_factor > 10

