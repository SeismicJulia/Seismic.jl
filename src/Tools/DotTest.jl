function DotTest(m_rand,d_rand,operators,parameters)
	# Dot product test for a vector of linear operators

	m_adj = LinearOperator(d_rand,operators,parameters,adj=true)
	d_fwd = LinearOperator(m_rand,operators,parameters,adj=false)
	inner1 = InnerProduct(d_rand,d_fwd)
	inner2 = InnerProduct(m_rand,m_adj)
	println("<d_rand,d_fwd> = ",inner1)
	println("<m_rand,m_adj> = ",inner2)
	return inner1,inner2
	
end

function DotTest(m_rand::ASCIIString,d_rand::ASCIIString,operators,parameters)
	# Dot product test for a vector of linear operators

	rand_string = string(round(Int,rand()*100000))
	m_adj = join(["tmp_DotTest_m_adj_",rand_string])
	d_fwd = join(["tmp_DotTest_d_adj_",rand_string])
	LinearOperator(m_adj,d_rand,operators,parameters,adj=true)
	LinearOperator(m_rand,d_fwd,operators,parameters,adj=false)
	inner1 = InnerProduct(d_rand,d_fwd)
	inner2 = InnerProduct(m_rand,m_adj)
	println("These values should be close to each other:")
	println("<d_rand,d_fwd> = ",inner1)
	println("<m_rand,m_adj> = ",inner2)
	SeisRemove(m_adj)
	SeisRemove(d_fwd)

end

function DotTest(m_rand::Array{ASCIIString,1},d_rand::Array{ASCIIString,1},operators,parameters)
	# Dot product test for a vector of linear operators
	# for 3 component data
	
	rand_string = string(round(Int,rand()*100000))
	m1_adj = join(["tmp_DotTest_m1_adj_",rand_string])
	m2_adj = join(["tmp_DotTest_m2_adj_",rand_string])
	d1_fwd = join(["tmp_DotTest_d1_adj_",rand_string])
	d2_fwd = join(["tmp_DotTest_d2_adj_",rand_string])
	m_adj = [m1_adj;m2_adj]
	d_fwd = [d1_fwd;d2_fwd]
	LinearOperator(m_adj,d_rand,operators,parameters,adj=true)
	LinearOperator(m_rand,d_fwd,operators,parameters,adj=false)
	inner1 = InnerProduct(d_rand,d_fwd)
	inner2 = InnerProduct(m_rand,m_adj)
	println("These values should be close to each other:")
	println("<d_rand,d_fwd> = ",inner1)
	println("<m_rand,m_adj> = ",inner2)
	SeisRemove(m_adj);
	SeisRemove(d_fwd);

end
