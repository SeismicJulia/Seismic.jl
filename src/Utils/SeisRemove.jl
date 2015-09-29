function SeisRemove(in::ASCIIString)

	rm(join([in ".seisd"]));rm(join([in ".seish"]));

end

function SeisRemove(in::Array{ASCIIString,1})

	for j = 1 : length(in)
		SeisRemove(in[j])
	end

end
