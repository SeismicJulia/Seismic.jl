function MatrixMultiplyOp(in,adj;matrix=1)

	if (adj)
		out = matrix'*in
	else
		out = matrix*in
	end

	return out
end

function MatrixMultiplyOp(in,h::Array{Header,1},e;adj=false,matrix=1)

	if (adj)
		out = zeros(Float32,size(matrix,2),size(in,2))
	else
		out = zeros(Float32,size(matrix,1),size(in,2))
	end
	for j = 1 : length(h)
		if (adj)
			out[:,j] = matrix'*in[:,j]
			h[j].n1 = size(matrix,2)
		else
			out[:,j] = matrix*in[:,j]
			h[j].n1 = size(matrix,1)
		end
	end
	return out,h
end

function MatrixMultiplyOp(m::AbstractString,d::AbstractString,adj;matrix=1)

	@compat parameters = Dict(:adj=>adj,:matrix=>matrix)
	if (adj==true)
		SeisProcess(d,m,[MatrixMultiplyOp],[parameters];key=["imx"])
		ext = ReadTextHeader(m)
		ext.n1 = size(matrix,2)
		Seismic.WriteTextHeader(m,ext,"native_float",4,join([m "@data@"]),join([m "@headers@"]))
	else
		SeisProcess(m,d,[MatrixMultiplyOp],[parameters];key=["imx"])
		ext = ReadTextHeader(d)
		ext.n1 = size(matrix,1)
		Seismic.WriteTextHeader(d,ext,"native_float",4,join([d "@data@"]),join([d "@headers@"]))
	end

end
