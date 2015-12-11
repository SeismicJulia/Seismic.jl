function IRLS(m0,d,param=Dict())
	# Non-quadratic regularization with Iteratively Reweighted Least Squares (IRLS).

	Niter_external = get(param,"Niter_external",2)
	Niter_internal = get(param,"Niter_internal",15)
	wm = get(param,"wm",1.)
	param["Niter"] = Niter_internal
	misfit = Float64[]
	m = copy(m0)
	v = copy(m0)
	param["wm"] = wm
	for iter_external = 1 : Niter_external
		if (length(Niter_internal) > 1)
			param["Niter"] = Niter_internal[iter_external]
		end
		v,misfit1 = ConjugateGradients(m0,d,param)
		append!(misfit,misfit1)
		m0 = copy(v)
		m = wm.*v
		wm = abs(m0./maximum(abs(m0[:])))
		param["wm"] = wm
	end

	return m, misfit

end
