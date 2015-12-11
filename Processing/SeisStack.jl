function stacktraces(d,h,param)
	normalize = get(param,"normalize",true)
	h_out = [h[1]]
	h_out[1].trid = size(d,2)
	if (normalize == true)
		val = sum(d,2)/size(d,2)
	else
		val = sum(d,2)
	end
	return val, h_out
end

function SeisStack(in::ASCIIString,out::ASCIIString,param=Dict())

	param["key"] = get(param,"key",["imx","imy"])
	param["group"] = "gather"
	param["f"] = [stacktraces]
	SeisProcess(in,out,param)

end
