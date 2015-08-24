function DotTest(m_rand,d_rand,op,param=Dict())
# Dot product test for a vector of linear operators

	m_adj = Seismic.adjoint_op(d_rand,op,param)
	d_fwd = Seismic.forward_op(m_rand,op,param)
	inner1 = InnerProduct(d_rand[:],d_fwd[:])
	inner2 = InnerProduct(m_rand[:],m_adj[:])
	println("These values should be close to each other:")
	println("dot(d_rand,d_fwd) = ",inner1)
	println("dot(m_rand,m_adj) = ",inner2)

end

function DotTest(m_rand::ASCIIString,d_rand::ASCIIString,op,param=Dict())
# Dot product test for a vector of linear operators

	Seismic.adjoint_op("m_adj",d_rand,op,param)
	Seismic.forward_op(m_rand,"d_fwd",op,param)
	inner1 = InnerProduct(d_rand,"d_fwd")
	inner2 = InnerProduct(m_rand,"m_adj")
	println("These values should be close to each other:")
	println("dot(d_rand,d_fwd) = ",inner1)
	println("dot(m_rand,m_adj) = ",inner2)

end
