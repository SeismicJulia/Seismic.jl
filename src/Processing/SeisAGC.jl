#simple AGC with a window specified in samples
#written by GRAM

function SeisAGC(in, param = Dict())

	nt = size(in,1)
	nx = size(in,2)
	out = zeros(size(in))

	wind = get(param,"wind",250)

	for i = 1 : nt, j = 1 : nx
		if i <= wind #slide on window
			out[i, j] = in[i, j] / abs(mean(in[1 : i, j]))
		else #full window
			out[i, j] = in[i, j] / abs(mean(in[i - wind : i, j]))
		end
	end

	return out

end
