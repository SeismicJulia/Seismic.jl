function ConjugateGradients(d,param=Dict())
# Non-quadratic regularization with CG-LS. The inner CG routine is taken from
# Algorithm 2 from Scales, 1987. Make sure linear operator passes the dot product.

	Niter_external = get(param,"Niter_external",1)
	Niter_internal = get(param,"Niter_internal",15)
	wd = get(param,"wd",ones(Float64,size(d))) # data weights
	operators = get(param,"operators",[fft_op])
	misfit = Array(Float64,1)
	global m
	for iter_external = 1 : Niter_external
		if (iter_external > 1)
			Pv = v.*P.*wm
			r = forward_op(Pv,operators,param)
			r = r.*wd
    			r = d - r
    		else
    			r = copy(d)
    		end	
    		g = adjoint_op(r,operators,param)
    		if (iter_external == 1)
    				wm = get(param,"wm",ones(Float64,size(g))) # model weights
    		end
    		g = g.*wm
    		if (iter_external == 1)
    			m = g
    			v = g
    			P = m*0. + 1.
    		else
    			g = g.*P
    		end
      		s = copy(g)
    		gamma = InnerProduct(g[:],g[:])
    		gamma_old = copy(gamma)
		for iter_internal = 1 : Niter_internal
			Ps = s.*P
			ss = forward_op(Ps,operators,param)
			ss = ss.*wd
			delta = InnerProduct(ss[:],ss[:])
      			alpha = gamma/(delta + 0.00000001)
			v = CGStep(v,s,1.,alpha)
			r = CGStep(r,ss,1.,-alpha)
			push!(misfit,InnerProduct(r[:],r[:]))
      			r = r.*wd
      			g = adjoint_op(r,operators,param)
      			g = g.*P.*wm
			gamma = InnerProduct(g[:],g[:])
			beta = gamma/(gamma_old + 0.00000001)
			gamma_old = copy(gamma)
      			s = CGStep(s,g,beta,1)
    		end
    		if (Niter_external > 1)
    			m = v.*P;
    			P = abs(m./maximum(abs(m)))
    		else
    			m = v;	
    		end
	end

	return m, misfit
end

function ConjugateGradients(m::ASCIIString,d::ASCIIString,misfit::ASCIIString,param=Dict())
# Non-quadratic regularization with CG-LS. The inner CG routine is taken from
# Algorithm 2 from Scales, 1987. Make sure linear operator passes the dot product.

	Niter_external = get(param,"Niter_external",1)
	Niter_internal = get(param,"Niter_internal",15)
	wd = get(param,"wd","tmp_cg_wd") # data weights
	operators = get(param,"operators",[fft_op])
	op = get(param,"op",fft_op)
	P = join(["tmp_CG_P_",string(int(rand()*100000))])
	v = join(["tmp_CG_v_",string(int(rand()*100000))])
	g = join(["tmp_CG_g_",string(int(rand()*100000))])
	s = join(["tmp_CG_s_",string(int(rand()*100000))])
	ss = join(["tmp_CG_ss_",string(int(rand()*100000))])
	r = join(["tmp_CG_r_",string(int(rand()*100000))])
	Pv = join(["tmp_CG_Pv_",string(int(rand()*100000))])
	Ps = join(["tmp_CG_Ps_",string(int(rand()*100000))])
	SeisCopy(d,r)
	fp = open(misfit,"w")
	write(fp,join(["began execution at: ",strftime(time()),"\n"]))
	write(fp,join([string(InnerProduct(r,r)),"\n"]))
	close(fp)

	for iter_external = 1 : Niter_external
		if (iter_external > 1)
			CGMult(Pv,P,v,1.,0.)
			CGMult(Pv,Pv,wm,1.,0.)
			forward_op(Pv,r,operators,param)
    			CGMult(r,r,wd,1.,0)
    			CGStep(r,d,-1.,1.,0.)
    		else
    			SeisCopy(d,r)
    		end	
    		adjoint_op(g,r,operators,param)
    		if (iter_external == 1)
			SeisCopy(g,m)
			SeisCopy(g,v)
    			if (!haskey(param,"wm"))
				wm = join(["tmp_CG_wm_",string(int(rand()*100000))])
				CGMult(wm,g,m,0.,1.)
			end
    		end
		CGMult(g,g,wm,1.,0.)
    		if (iter_external == 1)
			CGMult(m,m,wm,1.,0.)
			CGMult(v,v,wm,1.,0.)
    	    		CGMult(P,g,m,0.,1.)
    		else
    			CGMult(g,g,P,1.,0.)    
    		end
        	SeisCopy(g,s)
    		gamma = InnerProduct(g,g)
    		gamma_old = copy(gamma)
		for iter_internal = 1 : Niter_internal
			CGMult(Ps,s,P,1.,0.)
			forward_op(Ps,ss,operators,param)
			CGMult(ss,ss,wd,1.,0.)
			delta = InnerProduct(ss,ss)
      			alpha = gamma/(delta + 0.00000001)
			CGStep(v,s,1.,alpha)
			CGStep(r,ss,1.,-alpha)
			fp = open(misfit,"a")
			write(fp,join([string(InnerProduct(r,r)),"\n"]))
			close(fp)
      			CGMult(r,r,wd,1.,0.)
      			adjoint_op(g,r,operators,param)
      			CGMult(g,g,wm,1.,0.) 
      			CGMult(g,g,P,1.,0.) 
			gamma = InnerProduct(g,g)
			beta = gamma/(gamma_old + 0.00000001)
			gamma_old = copy(gamma)
      			CGStep(s,g,beta,1.)
    		end
    		if (Niter_external > 1)
    			CGMult(m,v,P,1.,0.)
    			CGSparseNorm(P,m)
    		else
    			SeisCopy(v,m)		
    		end
	end
	fp = open(misfit,"a")
	write(fp,join(["finished execution at: ",strftime(time()),"\n"]))
	close(fp)

end

function forward_op(m,operators,param)
	param["adj"] = false
	d = [];
	for iop = 1 : 1 : length(operators)
		op = operators[iop]
		d = op(m,param)
		m = copy(d)
	end
	return d
end

function adjoint_op(d,operators,param)
	param["adj"] = true
	m = [];
	for iop = length(operators) : -1 : 1
		op = operators[iop]
		m = op(d,param)
		d = copy(m)
	end
	return m
end

function forward_op(m::ASCIIString,d::ASCIIString,operators,param)
	param["adj"] = false
	tmp_filename_m = join(["tmp_CGFWD_m_",string(int(rand()*100000))])
	tmp_filename_d = join(["tmp_CGFWD_d_",string(int(rand()*100000))])
	SeisCopy(m,tmp_filename_m)
	for iop = 1 : 1 : length(operators)
		op = operators[iop]
		op(tmp_filename_m,tmp_filename_d,param)
		SeisCopy(tmp_filename_d,tmp_filename_m)
	end
	SeisCopy(tmp_filename_d,d)
	SeisRemove(tmp_filename_m)
	SeisRemove(tmp_filename_d)
end

function adjoint_op(m::ASCIIString,d::ASCIIString,operators,param)
	param["adj"] = true
	tmp_filename_m = join(["tmp_CGADJ_m_",string(int(rand()*100000))])
	tmp_filename_d = join(["tmp_CGADJ_d_",string(int(rand()*100000))])
	SeisCopy(d,tmp_filename_d)
	for iop = length(operators) : -1 : 1
		op = operators[iop]
		op(tmp_filename_m,tmp_filename_d,param)
		SeisCopy(tmp_filename_m,tmp_filename_d)
	end
	SeisCopy(tmp_filename_m,m)
	SeisRemove(tmp_filename_m)
	SeisRemove(tmp_filename_d)
end

