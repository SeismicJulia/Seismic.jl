function WeightingOp(in,adj;w=1)

	return in.*w

end

function WeightingOp(m::AbstractString,d::AbstractString,adj;w="NULL")

	if (adj==true)
		d1,h1,e1 = SeisRead(d)
		d2,h2,e2 = SeisRead(w)
		SeisWrite(m,d1[:,:].*d2[:,:],h1,e1)
	else
		d1,h1,e1 = SeisRead(m)
		d2,h2,e2 = SeisRead(w)
		SeisWrite(d,d1[:,:].*d2[:,:],h1,e1)
	end

end

function WeightingOp(m::Array{AbstractString,1},d::Array{AbstractString,1},adj;w="NULL")

	if !(isdefined(w,2))
		for j = 1 : length(m)
			WeightingOp(m[j],d[j],adj,w=w)
		end
	else
		for j = 1 : length(m)
			WeightingOp(m[j],d[j],adj,w=w[j])
		end
	end
end
