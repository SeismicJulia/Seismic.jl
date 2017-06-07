function SeisSmooth1(A;L=7)
	a = isodd(L) ? 0 : 1
	kernel = ones(eltype(A),L)/L;
	A = conv(A,kernel)
	A = A[floor(Int,L/2) + 1:end-floor(Int,L/2) + a]
	return A
end

function SeisSmooth1(A,h::Array{Header,1};L=7)
	for itrace = 1 : size(A[:,:],2)
		A[:,itrace] = SeisSmooth1(A[:,itrace],L=L)
	end
	return A,h
end

function SeisSmooth1(in::AbstractString,out::AbstractString;L=7)
	@compat parameters = Dict(:L=>L)
	SeisProcess(in,out,[SeisSmooth1],[parameters];group="some",ntrace=100000)
end
