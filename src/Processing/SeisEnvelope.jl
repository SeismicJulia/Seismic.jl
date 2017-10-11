function SeisEnvelope(d)

	D = fft(d,1)
	D[1:Int(floor(size(d,1)/2)),:] = 0.0
	return 2*abs.(ifft(D,1))

end

function SeisEnvelope(A,h::Array{Header,1})
	for itrace = 1 : size(A[:,:],2)
		A[:,itrace] = SeisEnvelope(A[:,itrace])
	end
	return A,h
end

function SeisEnvelope(in::AbstractString,out::AbstractString)
	@compat parameters = Dict()
	SeisProcess(in,out,[SeisEnvelope],[parameters];group="some",ntrace=100000)
end
