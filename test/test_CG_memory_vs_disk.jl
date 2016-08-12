using Seismic
using Base.Test

# test of Conjugate Gradients
nd = 100
nm = 300
L = 10*randn(nd,nm)
m = 10.0*rand(nm)
d = zeros(nd)
d = L*m
d = d + 2.0*randn(nd)

println(norm(d))
println(norm(m))

m1,cost = ConjugateGradients(d,[MatrixMultiplyOp],[Dict(:matrix=>L)],Niter=20,mu=0.2)

ext = Seismic.Extent(size(d,1),size(d,2),1,1,1,
0,0,0,1,1,
1,1,1,1,1,
"","","","","",
"","","","","",
"")
h = Header[]
for ix = 1:size(d,2)
		push!(h,Seismic.InitSeisHeader())
		h[ix] = Seismic.InitSeisHeader()
		h[ix].tracenum = ix
		h[ix].n1 = size(d,1)
		h[ix].d1 = 1
end

SeisWrite("tmp_d.seis",d,h,ext)

ConjugateGradients("tmp_m.seis","tmp_d.seis",[MatrixMultiplyOp],[Dict(:matrix=>L)],"tmp_cost.txt",Niter=20,mu=0.2)
m2,ext = SeisRead("tmp_m.seis")

# test that quality factor between disk and memory based CG 
# is greater than 50 Decibels
quality_factor = 10*log10(norm(m1[:],2)/norm(m2[:]-m1[:],2))
println("Quality factor = ",quality_factor)
@test quality_factor > 10.
