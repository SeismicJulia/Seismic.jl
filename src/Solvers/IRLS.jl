function IRLS(d,operators,parameters;Niter_external=3,Niter_internal=10,mu=0)
	# Non-quadratic regularization with Iteratively Reweighted Least Squares (IRLS).

	cost = Float64[]
	weights = ones(Float64,size(d))
	parameters[end][:w] = weights
	m = []
	for iter_external = 1 : Niter_external
		v,cost1 = ConjugateGradients(d,operators,parameters,Niter=Niter_internal,mu=mu)
		append!(cost,cost1)
		m = v.*weights
		weights = abs.(m./maximum(abs.(m[:])))
		parameters[end][:w] = weights
	end

	return m, cost

end
