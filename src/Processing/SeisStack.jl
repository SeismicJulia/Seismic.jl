function stacktraces(d,h,param)
	h_out = [h[1]]
	h_out[1].trid = size(d,2)
	return sum(d,2)/size(d,2), h_out
end

function SeisStack(in::ASCIIString,out::ASCIIString,param)

	key = get(param,"key",["imx","imy"])
	SeisProcess(in,out,[stacktraces],param,"gather",key)

end
