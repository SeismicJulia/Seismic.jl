function stacktraces(d,h;normalize=true)
	h_out = [h[1]]
	h_out[1].trid = size(d[:,:],2)
	if (normalize == true)
		val = sum(d[:,:],2)/size(d[:,:],2)
	else
		val = sum(d[:,:],2)
	end
	return val, h_out
end

function SeisStack(in::ASCIIString,out::ASCIIString;key=["imx" "imy"],normalize=true)

	@compat parameters = Dict(:normalize=>true)
	SeisProcess(in,out,[stacktraces],[parameters],key=key)
	DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
	filename_data_out = join([DATAPATH out "@data@"])
	filename_headers_out = join([DATAPATH out "@headers@"])
	extent = ReadTextHeader(in)

	stream_h = open(filename_headers_out)
	@compat extent.n2 = round(Int,filesize(stream_h)/(4*length(fieldnames(Header))))
	close(stream_h)
	extent.d2 = 1.; extent.o2 = 0.; extent.label2 = ""; extent.unit2 = "";
	extent.n3 = 1; extent.d3 = 1.; extent.o3 = 0.; extent.label3 = ""; extent.unit3 = "";
	extent.n4 = 1; extent.d4 = 1.; extent.o4 = 0.; extent.label4 = ""; extent.unit4 = "";
	extent.n5 = 1; extent.d5 = 1.; extent.o5 = 0.; extent.label5 = ""; extent.unit5 = "";
	WriteTextHeader(out,extent,"native_float",4,filename_data_out,filename_headers_out)

end
