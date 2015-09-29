function DotTest(m_rand,d_rand,op,param=Dict())
	# Dot product test for a vector of linear operators

	m_adj = Seismic.adjoint_op(d_rand,op,param)
	d_fwd = Seismic.forward_op(m_rand,op,param)
	inner1 = InnerProduct(d_rand[:],d_fwd[:])
	inner2 = InnerProduct(m_rand[:],m_adj[:])
	println("These values should be close to each other:")
	println("<d_rand,d_fwd> = ",inner1)
	println("<m_rand,m_adj> = ",inner2)

end

function DotTest(m_rand::ASCIIString,d_rand::ASCIIString,op,param=Dict())
	# Dot product test for a vector of linear operators

	m_adj = join(["tmp_DotTest_m_adj_",string(int(rand()*100000))])
	d_fwd = join(["tmp_DotTest_d_adj_",string(int(rand()*100000))])
	Seismic.adjoint_op(m_adj,d_rand,op,param)
	Seismic.forward_op(m_rand,d_fwd,op,param)
	inner1 = InnerProduct(d_rand,d_fwd)
	inner2 = InnerProduct(m_rand,m_adj)
	println("These values should be close to each other:")
	println("<d_rand,d_fwd> = ",inner1)
	println("<m_rand,m_adj> = ",inner2)
	SeisRemove(m_adj)
	SeisRemove(d_fwd)

end

function DotTest(m1_rand::ASCIIString,m2_rand::ASCIIString,m3_rand::ASCIIString,d1_rand::ASCIIString,d2_rand::ASCIIString,d3_rand::ASCIIString,op,param=Dict())
	# Dot product test for a vector of linear operators
	# for 3 component data
	
	m1_adj = join(["tmp_DotTest_m1_adj_",string(int(rand()*100000))])
	m2_adj = join(["tmp_DotTest_m2_adj_",string(int(rand()*100000))])
	m3_adj = join(["tmp_DotTest_m3_adj_",string(int(rand()*100000))])
	d1_fwd = join(["tmp_DotTest_d1_adj_",string(int(rand()*100000))])
	d2_fwd = join(["tmp_DotTest_d2_adj_",string(int(rand()*100000))])
	d3_fwd = join(["tmp_DotTest_d3_adj_",string(int(rand()*100000))])
	Seismic.adjoint_op(m1_adj,m2_adj,m3_adj,d1_rand,d2_rand,d3_rand,op,param)
	Seismic.forward_op(m1_rand,m2_rand,m3_rand,d1_fwd,d2_fwd,d3_fwd,op,param)
	inner1 = InnerProduct(d1_rand,d1_fwd) + InnerProduct(d2_rand,d2_fwd) + InnerProduct(d3_rand,d3_fwd)
	inner2 = InnerProduct(m1_rand,m1_adj) + InnerProduct(m2_rand,m2_adj) + InnerProduct(m3_rand,m3_adj)
	println("These values should be close to each other:")
	println("<d_rand,d_fwd> = ",inner1)
	println("<m_rand,m_adj> = ",inner2)
	SeisRemove(m1_adj);SeisRemove(m2_adj);SeisRemove(m3_adj)
	SeisRemove(d1_fwd);SeisRemove(d2_fwd);SeisRemove(d3_fwd)

end
