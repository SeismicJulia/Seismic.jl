function DotTest(m_rand,d_rand,op,param=Dict())
# Dot product test for a linear operator

	param["adj"] = true
	m_adj = op(d_rand,param)
	param["adj"] = false
	d_fwd = op(m_rand,param)
	inner1 = InnerProduct(d_rand[:],d_fwd[:])
	inner2 = InnerProduct(m_rand[:],m_adj[:])
	println("These values should be close to each other:")
	println("dot(d_rand,d_fwd) = ",inner1)
	println("dot(m_rand,m_adj) = ",inner2)

end

function DotTest(m_rand::ASCIIString,d_rand::ASCIIString,op,param=Dict())
# Dot product test for a linear operator

	param["adj"] = true
	op("m_adj",d_rand,param)
	param["adj"] = false
	op(m_rand,"d_fwd",param)
	inner1 = InnerProduct(d_rand,"d_fwd")
	inner2 = InnerProduct(m_rand,"m_adj")
	println("These values should be close to each other:")
	println("dot(d_rand,d_fwd) = ",inner1)
	println("dot(m_rand,m_adj) = ",inner2)

end
