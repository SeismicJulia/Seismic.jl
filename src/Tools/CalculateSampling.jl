function CalculateSampling(in)

	cutoff = 1e-10
	itrace = 1
	
	n=size(in)
	in=reshape(in,n[1],:)
        wd = zeros(Float32,size(in))	
	n2=size(in,2)
	for itrace = 1 : n2
		a = sum(in[:,itrace].*in[:,itrace])
		if (a > cutoff)
			wd[:,itrace] = 1.
		end
	end
	wd=reshape(wd,n)
	return wd;
end 

function CalculateSampling(in,h;cutoff=1e-10)

	itrace = 1
        n=size(in)
	in=reshape(in,n[1],:)
	wd = zeros(Float32,size(in))
	n2=size(in,2)
	for itrace = 1 : n2
		a = sqrt(sum(in[:,itrace].^2))
		if (a > cutoff)
			wd[:,itrace] = 1.
		end
	end
	wd=reshape(wd,n)
	return wd,h;
end






function CalculateSampling(in::AbstractString,wd::AbstractString;cutoff=1e-10)
	# calculate sampling operator (1's for live traces, 0's for missing traces)

	@compat parameters = Dict(:cutoff=>cutoff)
	SeisProcess(in,wd,[CalculateSampling],[parameters])

end
