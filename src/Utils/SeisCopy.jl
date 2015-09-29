function SeisCopy(in::ASCIIString,out::ASCIIString)

	run(`cp ${join([in ".seisd"])} ${join([out ".seisd"])}`);run(`cp ${join([in ".seish"])} ${join([out ".seish"])}`);

end

function SeisCopy(in::Array{ASCIIString,1},out::Array{ASCIIString,1})

	for j = 1 : length(in)
		SeisCopy(in[j],out[j])
	end
	
end
