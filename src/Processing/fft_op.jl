function fft_op(in,param=Dict())
	adj = get(param,"adj",false)
	if (adj)
		out = fft(in)/sqrt(length(in[:]))
	else
		out = bfft(in)/sqrt(length(in[:]))	
	end
	
	return out
end
