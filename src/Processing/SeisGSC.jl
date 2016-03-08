#geometrical spreading correction by sample * time^power
#written by GRAM

function SeisGSC(in,param=Dict())

	nt = size(in,1)
	nx = size(in,2)
	out = zeros(size(in))

	tpow = get(param,"tpow",2)

	for i = 1 : nt, j = 1 : nx
		out[i,j] = in[i,j] * i^tpow
	end

	return out

end
