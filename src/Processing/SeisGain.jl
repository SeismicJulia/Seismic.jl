#simple data scale by multiplicative scalar
#written by GRAM

function SeisGain(in,param=Dict())

	nt = size(in,1)
	nx = size(in,2)
	out = zeros(size(in))

	scale = get(param,"scale",2)

	out = in[:] * scale

	return out

end
