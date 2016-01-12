using Seismic
using Base.Test

# test of Conjugate Gradients
nd = 1000
nm = 300
L = rand(nd,nm)
m = rand(nm);
d = L*m;
d = d + 0.2*rand(nd);

m0 = zeros(nm)
m1,cost = ConjugateGradients(m0,d,[[MatrixMultiplyOp]],[[:matrix_multiply=>0.1]],Niter=nm)

SeisWrite("d",d)

SeisWrite("m0",0.*m)

ConjugateGradients("m2","m0","d","misfit.txt",param)
m2,h = SeisRead("m2")

# test that quality factor between disk and memory based CG 
# is greater than 50 Decibels
quality_factor = 10*log10(norm(m1[:],2)/norm(m2[:]-m1[:],2))
println("Quality factor = ",quality_factor)
@test quality_factor > 10.
