using Seismic
using Base.Test

# test of Conjugate Gradients
nd = 1000
nm = 300
L = randn(nd,nm)
m = randn(nm)
d = L*m
d = d + 0.2*randn(nd)

# closed form solution
mu = 0.2
m1 = (L'*L + mu*eye(nm))\(L'*d)

m0 = zeros(nm)
m2,cost = ConjugateGradients(d,[MatrixMultiplyOp],[Dict(:matrix=>L)],Niter=nm,mu=mu)

# test that quality factor between CG and closed form solution 
# is greater than 50 Decibels
quality_factor = 10*log10(norm(m1[:],2)/norm(m2[:]-m1[:],2))
println("Quality factor = ",quality_factor)
@test quality_factor > 10.
