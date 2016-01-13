function WeightingOp(in,adj;w=1)
	
	return in.*w	

end

function WeightingOp(m::ASCIIString,d::ASCIIString,adj;w="NULL")

	if (adj==true)
		d1,h1,e1 = SeisRead(d)
		d2,h2,e2 = SeisRead(w)
		SeisWrite(m,d1.*d2,h1,e1)
	else
		d1,h1,e1 = SeisRead(m)
		d2,h2,e2 = SeisRead(w)
		SeisWrite(d,d1.*d2,h1,e1)
	end

end

function ApplyDataWeights(m::Array{ASCIIString,1},d::Array{ASCIIString,1},adj;w="NULL")

	for j = 1 : length(m)
		ApplyDataWeights(m[j],d[j],adj,w=w)
	end     

end                                                                                                                     

