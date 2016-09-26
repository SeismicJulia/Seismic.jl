using Seismic
using Base.Test

# test that FFTOp passes the dot product test
m_rand = randn(20,20,20,20) + 1im*randn(20,20,20,20);
d_rand = randn(20,20,20,20) + 1im*randn(20,20,20,20);
T = rand(20,20,20,20); a = find(T.<=0.5); b = find(T.>0.5); T[a] = 0.; T[b] = 1.;
operators = [WeightingOp,FFTOp,WeightingOp]
parameters = [Dict(:w=>T),Dict(:normalize=>true),Dict(:w=>ones(T))]
a,b = Seismic.DotTest(m_rand,d_rand,operators,parameters)
@test abs((a-b)/a) < 0.01

