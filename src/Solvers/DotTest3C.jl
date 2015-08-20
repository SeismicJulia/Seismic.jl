function DotTest3C(op,param,d_rand::Array{Float64,2},m_rand::Array{Float64,2})
# Dot product test for a linear operator for 3 component data (and 2 component model space)

	param["adj"] = "y"
	m_adj = op(d_rand,param)
	param["adj"] = "n"
	d_fwd = op(m_rand,param)
	inner1 = sum(conj(d_rand[:]).*d_fwd[:])
	inner2 = sum(conj(m_rand[:]).*m_adj[:])
	println("These values should be close to each other:")
	println("dot(d_rand,d_fwd) = ",inner1)
	println("dot(m_rand,m_adj) = ",inner2)

end

function DotTest3C(op,param,ux_rand::ASCIIString,uy_rand::ASCIIString,uz_rand::ASCIIString,mpp_rand::ASCIIString,mps_rand::ASCIIString)
# Dot product test for a linear operator for 3 component data (and 2 component model space)

	param["adj"] = "y"
	op(ux_rand,uy_rand,uz_rand,"mpp_adj","mps_adj",param)
	param["adj"] = "n"
	op("ux_fwd","uy_fwd","uz_fwd",mpp_rand,mps_rand,param)
	ux1,h,s = SeisRead(ux_rand)
	ux2,h,s = SeisRead("ux_fwd")
	uy1,h,s = SeisRead(uy_rand)
	uy2,h,s = SeisRead("uy_fwd")
	uz1,h,s = SeisRead(uz_rand)
	uz2,h,s = SeisRead("uz_fwd")
	mpp1,h,s = SeisRead(mpp_rand)
	mpp2,h,s = SeisRead("mpp_adj")
	mps1,h,s = SeisRead(mps_rand)
	mps2,h,s = SeisRead("mps_adj")
	println(sum(ux1[:].*ux1[:]))
	println(sum(ux2[:].*ux2[:]))
	println(sum(uy1[:].*uy1[:]))
	println(sum(uy2[:].*uy2[:]))
	println(sum(uz1[:].*uz1[:]))
	println(sum(uz2[:].*uz2[:]))
	inner1 = sum(conj(ux1[:]).*ux2[:]) + sum(conj(uy1[:]).*uy2[:]) + sum(conj(uz1[:]).*uz2[:])
	inner2 = sum(conj(mpp1[:]).*mpp2[:]) + sum(conj(mps1[:]).*mps2[:])
	println("These values should be close to each other:")
	println("dot(ux_rand,ux_fwd) + dot(uy_rand,uy_fwd) + dot(uz_rand,uz_fwd) = ",inner1)
	println("dot(mpp_rand,mpp_adj) + dot(mps_rand,mps_adj) = ",inner2)

end
