using Seismic
using Base.Test

# test of Conjugate Gradients
nd = 1000
nm = 300
L = rand(nd,nm)
m = rand(nm);
d = L*m;
d = d + 0.2*rand(nd);

# closed form solution
mu = 0.2;
m1 = (L'*L + mu*eye(nm))\(L'*d);

m0 = zeros(nm)
param = ["Niter"=>nm,"operators"=>[MatrixMultiply],"matrix"=>L,"mu"=>mu]
m2,cost = ConjugateGradients(m0,d,param)

# test that quality factor between CG and closed form solution 
# is greater than 50 Decibels
quality_factor = 10*log10(norm(m1[:],2)/norm(m2[:]-m1[:],2))
println("Quality factor = ",quality_factor)
@test quality_factor > 50.

