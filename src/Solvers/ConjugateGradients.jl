function ConjugateGradients(m0,d,param::Dict{Any,Any})
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	Niter = get(param,"Niter",10)
	mu = get(param,"mu",0.)
	mu = sqrt(mu)
	operators = get(param,"operators",[])
	param_op = copy(param)
	cost = Float64[]

	m = copy(m0)
	r = forward_op(m0,operators,param_op)
	r = d - r
	rr = 0. - mu*m0;
	push!(cost,InnerProduct(r[:],r[:]) + InnerProduct(rr[:],rr[:]))
	g = adjoint_op(r,operators,param_op)
	g = g + mu.*rr
	s = copy(g)
	gamma_old = InnerProduct(g[:],g[:])
	for iter = 1 : Niter
		t = forward_op(s,operators,param_op)
		tt = mu.*s
		delta = InnerProduct(t[:],t[:]) + InnerProduct(tt[:],tt[:])
		alpha = gamma_old/(delta + 1.e-20)
		m = m + alpha*s
		r = r - alpha*t
		rr = rr - alpha*tt
		push!(cost,InnerProduct(r[:],r[:]) + InnerProduct(rr[:],rr[:]))
		g = adjoint_op(r,operators,param_op)
		g = g + mu*rr
		gamma = InnerProduct(g[:],g[:])
		beta = gamma/(gamma_old + 1.e-20)
		gamma_old = copy(gamma)
		s = beta*s + g
	end

	return m, cost
end

function ConjugateGradients(m::ASCIIString,m0::ASCIIString,d::ASCIIString,cost_file::ASCIIString,param::Dict{Any,Any})
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	Niter = get(param,"Niter",10)
	mu = get(param,"mu",0.)
	mu = sqrt(mu)
	operators = get(param,"operators",[])
	param_op = copy(param)
	cost = Float64[]
	rand_string = string(int(rand()*100000))
	g = join(["tmp_CG_g_",rand_string])
	s = join(["tmp_CG_s_",rand_string])
	rr = join(["tmp_CG_rr_",rand_string])
	tt = join(["tmp_CG_tt_",rand_string])
	r = join(["tmp_CG_r_",rand_string])
	t = join(["tmp_CG_t_",rand_string])
	SeisCopy(m0,m)
	forward_op(m0,r,operators,param_op)
	CGStep(r,d,{"a"=>-1.,"b"=>1.})
	SeisCopy(m0,rr)	
	CGStep(rr,m0,{"a"=>0.,"b"=>-mu})
	push!(cost,InnerProduct(r,r) + InnerProduct(rr,rr))
	fp = open(cost_file,"w")
	write(fp,join(["began execution at: ",strftime(time()),"\n"]))
	write(fp,join([string(cost[1]),"\n"]))
	close(fp)
	adjoint_op(g,r,operators,param_op)
	CGStep(g,rr,{"a"=>1.,"b"=>mu})
	SeisCopy(g,s)
	# initialize tt with zeros
	SeisCopy(s,tt)
	CGStep(tt,s,{"a"=>0.,"b"=>0.})
	gamma_old = InnerProduct(g,g)
	for iter = 1 : Niter	
		forward_op(s,t,operators,param_op)
		CGStep(tt,s,{"a"=>0.,"b"=>mu})
		delta = InnerProduct(t,t) + InnerProduct(tt,tt)
		alpha = gamma_old/(delta + 1.e-20)
		CGStep(m,s,{"a"=>1.,"b"=>alpha})
		CGStep(r,t,{"a"=>1.,"b"=>-alpha})
		CGStep(rr,tt,{"a"=>1.,"b"=>-alpha})
		push!(cost,InnerProduct(r,r) + InnerProduct(rr,rr))
		fp = open(cost_file,"a")
		write(fp,join([string(cost[iter+1]),"\n"]))
		close(fp)
		adjoint_op(g,r,operators,param_op)
		CGStep(g,rr,{"a"=>1.,"b"=>mu})
		gamma = InnerProduct(g,g)
		println("gamma=",gamma)
		beta = gamma/(gamma_old + 1.e-20)
		gamma_old = copy(gamma)
		CGStep(s,g,{"a"=>beta,"b"=>1.})
	end
	SeisRemove(g);
	SeisRemove(s);
	SeisRemove(rr);
	SeisRemove(tt)
	SeisRemove(r);
	SeisRemove(t);

end

function ConjugateGradients(m::Array{ASCIIString,1},m0::Array{ASCIIString,1},d::Array{ASCIIString,1},cost_file::ASCIIString,param::Dict{Any,Any})
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	Niter = get(param,"Niter",10)
	mu = get(param,"mu",[0. 0. 0.])
	mu = sqrt(mu)
	operators = get(param,"operators",[])
	param_op = copy(param)
	cost = Float64[]
	rand_string = string(int(rand()*100000))
	g = [join(["tmp_CG_g1_",rand_string]);join(["tmp_CG_g2_",rand_string]);join(["tmp_CG_g3_",rand_string])]
	s = [join(["tmp_CG_s1_",rand_string]);join(["tmp_CG_s2_",rand_string]);join(["tmp_CG_s3_",rand_string])]
	rr = [join(["tmp_CG_rr1_",rand_string]);join(["tmp_CG_rr2_",rand_string]);join(["tmp_CG_rr3_",rand_string])]
	tt = [join(["tmp_CG_tt1_",rand_string]);join(["tmp_CG_tt2_",rand_string]);join(["tmp_CG_tt3_",rand_string])]
	r = [join(["tmp_CG_r1_",rand_string]);join(["tmp_CG_r2_",rand_string]);join(["tmp_CG_r3_",rand_string])]
	t = [join(["tmp_CG_t1_",rand_string]);join(["tmp_CG_t2_",rand_string]);join(["tmp_CG_t3_",rand_string])]
	SeisCopy(m0,m)
	forward_op(m0,r,operators,param_op)
	CGStep(r,d,{"a"=>[-1.;-1.;-1.],"b"=>[1.;1.;1.]});
	SeisCopy(m0,rr);	
	CGStep(rr,m0,{"a"=>[0.;0.;0.],"b"=>-mu})
	push!(cost,InnerProduct(r,r) + InnerProduct(rr,rr))
	fp = open(cost_file,"w")
	write(fp,join(["began execution at: ",strftime(time()),"\n"]))
	write(fp,join([string(cost[1]),"\n"]))
	close(fp)
	adjoint_op(g,r,operators,param_op)
	CGStep(g,rr,{"a"=>[1.;1.;1.],"b"=>mu})
	SeisCopy(g,s)
	# initialize tt with zeros
	SeisCopy(s,tt)
	CGStep(tt,s,{"a"=>[0.;0.;0.],"b"=>[0.;0.;0.]});
	gamma_old = InnerProduct(g,g)	
	for iter = 1 : Niter	
		forward_op(s,t,operators,param_op)
		CGStep(tt,s,{"a"=>[0.;0.;0.],"b"=>mu});
		delta = InnerProduct(t,t) + InnerProduct(tt,tt)
		alpha = gamma_old/(delta + 1.e-20)
		CGStep(m,s,{"a"=>[1.;1.;1.],"b"=>[alpha;alpha;alpha]})
		CGStep(r,t,{"a"=>[1.;1.;1.],"b"=>-[alpha;alpha;alpha]})
		CGStep(rr,tt,{"a"=>[1.;1.;1.],"b"=>-[alpha;alpha;alpha]})
		push!(cost,InnerProduct(r,r) + InnerProduct(rr,rr))
		fp = open(cost_file,"a")
		write(fp,join([string(cost[iter+1]),"\n"]))
		close(fp)
		adjoint_op(g,r,operators,param_op)
		CGStep(g,rr,{"a"=>[1.;1.;1.],"b"=>mu})
		gamma = InnerProduct(g,g)
		beta = gamma/(gamma_old + 1.e-20)
		gamma_old = copy(gamma)
		CGStep(s,g,{"a"=>[beta;beta;beta],"b"=>[1.;1.;1.]})
	end
	SeisRemove(g);
	SeisRemove(s);
	SeisRemove(rr);
	SeisRemove(tt);
	SeisRemove(r);
	SeisRemove(t);

end

function forward_op(m,operators,param)
	param["adj"] = false
	d = [];
	if length(operators) > 0
		for iop = length(operators) : -1 : 1
			op = operators[iop]
			d = op(m,param)
			m = copy(d)
		end
	else
		d = copy(m)
	end
	return d
end

function adjoint_op(d,operators,param)
	param["adj"] = true
	m = [];
	if length(operators) > 0
		for iop = 1 : 1 : length(operators)
			op = operators[iop]
			m = op(d,param)
			d = copy(m)
		end
	else
		m = copy(d)
	end
	return m
end

function forward_op(m::ASCIIString,d::ASCIIString,operators,param)
	param["adj"] = false
	rand_string = string(int(rand()*100000))
	tmp_m = join(["tmp_CGADJ_m_",rand_string])
	tmp_d = join(["tmp_CGADJ_d_",rand_string])
	SeisCopy(m,tmp_m)
	for iop = length(operators) : -1 : 1
		op = operators[iop]
		op(tmp_m,tmp_d,param)
		SeisCopy(tmp_d,tmp_m)
	end
	SeisCopy(tmp_d,d)
	SeisRemove(tmp_m)
	SeisRemove(tmp_d)
end

function adjoint_op(m::ASCIIString,d::ASCIIString,operators,param)
	param["adj"] = true
	rand_string = string(int(rand()*100000))
	tmp_m = join(["tmp_CGADJ_m_",rand_string])
	tmp_d = join(["tmp_CGADJ_d_",rand_string])
	SeisCopy(d,tmp_d)
	for iop = 1 : 1 : length(operators)
		op = operators[iop]
		op(tmp_m,tmp_d,param)
		SeisCopy(tmp_m,tmp_d)
	end
	SeisCopy(tmp_m,m)
	SeisRemove(tmp_m)
	SeisRemove(tmp_d)
end

function forward_op(m::Array{ASCIIString,1},d::Array{ASCIIString,1},operators,param)
	param["adj"] = false
	rand_string = string(int(rand()*100000))
	tmp_m = [join(["tmp_CGADJ_m1_",rand_string]);join(["tmp_CGADJ_m2_",rand_string]);join(["tmp_CGADJ_m3_",rand_string])]
	tmp_d = [join(["tmp_CGADJ_d1_",rand_string]);join(["tmp_CGADJ_d2_",rand_string]);join(["tmp_CGADJ_d3_",rand_string])]
	SeisCopy(m,tmp_m)
	for iop = length(operators) : -1 : 1
		op = operators[iop]
		if (length(methods(op,(Array,Array,Dict{Any,Any}))) > 0)
			op(tmp_m,tmp_d,param)
			SeisCopy(tmp_d,tmp_m)
		else
			for j = 1 : length(m)
				op(tmp_m[j],tmp_d[j],param)
			end	
			SeisCopy(tmp_d,tmp_m)
		end
	end	
	SeisCopy(tmp_d,d)
	SeisRemove(tmp_m)
	SeisRemove(tmp_d)
end

function adjoint_op(m::Array{ASCIIString,1},d::Array{ASCIIString,1},operators,param)
	param["adj"] = true
	rand_string = string(int(rand()*100000))
	tmp_m = [join(["tmp_CGADJ_m1_",rand_string]);join(["tmp_CGADJ_m2_",rand_string]);join(["tmp_CGADJ_m3_",rand_string])]
	tmp_d = [join(["tmp_CGADJ_d1_",rand_string]);join(["tmp_CGADJ_d2_",rand_string]);join(["tmp_CGADJ_d3_",rand_string])]
	SeisCopy(d,tmp_d)
	for iop = 1 : 1 : length(operators)
		op = operators[iop]
		if (length(methods(op,(Array,Array,Dict{Any,Any}))) > 0)
			op(tmp_m,tmp_d,param)
			SeisCopy(tmp_m,tmp_d)
		else
			for j = 1 : length(m)
				op(tmp_m[j],tmp_d[j],param)
			end	
			SeisCopy(tmp_m,tmp_d)
		end
	end
	SeisCopy(tmp_m,m)
	SeisRemove(tmp_m)
	SeisRemove(tmp_d)
end
