function ConjugateGradients(d,operators,parameters;Niter=10,mu=0)
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	mu = sqrt(mu)
	cost = Float64[]
	r = copy(d)
	push!(cost,InnerProduct(r,r))
	g = LinearOperator(r,operators,parameters,adj=true)
	m = zeros(g)
	rr = zeros(g)
	s = copy(g)
	gamma_old = InnerProduct(g,g)
	for iter = 1 : Niter
		t = LinearOperator(s,operators,parameters,adj=false)
		tt = mu.*s
		delta = InnerProduct(t,t) + InnerProduct(tt,tt)
		alpha = gamma_old/(delta + 1.e-20)
		m = m + alpha*s
		r = r - alpha*t
		rr = rr - alpha*tt
		push!(cost,InnerProduct(r,r) + InnerProduct(rr,rr))
		g = LinearOperator(r,operators,parameters,adj=true)
		g = g + mu*rr
		gamma = InnerProduct(g,g)
		beta = gamma/(gamma_old + 1.e-20)
		gamma_old = copy(gamma)
		s = beta*s + g
	end

	return m, cost
end

function ConjugateGradients(m::ASCIIString,d::ASCIIString,operators,parameters,cost_file::ASCIIString;Niter=10,mu=0)
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	mu = sqrt(mu)
	cost = Float64[]
	rand_string = string(int(rand()*100000))
	g = join(["tmp_CG_g_",rand_string])
	s = join(["tmp_CG_s_",rand_string])
	rr = join(["tmp_CG_rr_",rand_string])
	tt = join(["tmp_CG_tt_",rand_string])
	r = join(["tmp_CG_r_",rand_string])
	t = join(["tmp_CG_t_",rand_string])
    SeisCopy(d,r)
	push!(cost,InnerProduct(r,r))
	fp = open(cost_file,"w")
	write(fp,join(["began execution at: ",strftime(time()),"\n"]))
	write(fp,join([string(cost[1]),"\n"]))
	close(fp)
    LinearOperator(g,r,operators,parameters,adj=true)
	SeisCopy(g,s)
	SeisCopy(s,tt)
	CGStep(tt,s,a=0.,b=0.)
    SeisCopy(tt,m)
	
	gamma_old = InnerProduct(g,g)
	for iter = 1 : Niter	
        LinearOperator(s,t,operators,parameters,adj=false)
		CGStep(tt,s,a=0.,b=mu)
		delta = InnerProduct(t,t) + InnerProduct(tt,tt)
		alpha = gamma_old/(delta + 1.e-20)
		CGStep(m,s,a=1.,b=alpha)  
		CGStep(r,t,a=1.,b=-alpha)
		CGStep(rr,tt,a=1.,b=-alpha)
		push!(cost,InnerProduct(r,r) + InnerProduct(rr,rr))
		fp = open(cost_file,"a")
		write(fp,join([string(cost[iter+1]),"\n"]))
		close(fp)
        LinearOperator(g,r,operators,parameters,adj=true)
		CGStep(g,rr,a=1.,b=mu)
		gamma = InnerProduct(g,g)
		println("gamma=",gamma)
		beta = gamma/(gamma_old + 1.e-20)
		gamma_old = copy(gamma)
		CGStep(s,g,a=beta,b=1.)
	end
	SeisRemove(g);
	SeisRemove(s);
	SeisRemove(rr);
	SeisRemove(tt)
	SeisRemove(r);
	SeisRemove(t);

end

function ConjugateGradients(m::Array{ASCIIString,1},d::Array{ASCIIString,1},operators,parameters,cost_file::ASCIIString;Niter=10,mu=[0 0 0])
    # Conjugate Gradients following Algorithm 2 from Scales, 1987. 
    # The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

    mu = sqrt(mu)
    cost = Float64[]
    rand_string = string(int(rand()*100000))
    g = [join(["tmp_CG_g1_",rand_string]);join(["tmp_CG_g2_",rand_string]);join(["tmp_CG_g3_",rand_string])]
    s = [join(["tmp_CG_s1_",rand_string]);join(["tmp_CG_s2_",rand_string]);join(["tmp_CG_s3_",rand_string])]
    rr = [join(["tmp_CG_rr1_",rand_string]);join(["tmp_CG_rr2_",rand_string]);join(["tmp_CG_rr3_",rand_string])]
    tt = [join(["tmp_CG_tt1_",rand_string]);join(["tmp_CG_tt2_",rand_string]);join(["tmp_CG_tt3_",rand_string])]
    r = [join(["tmp_CG_r1_",rand_string]);join(["tmp_CG_r2_",rand_string]);join(["tmp_CG_r3_",rand_string])]
    t = [join(["tmp_CG_t1_",rand_string]);join(["tmp_CG_t2_",rand_string]);join(["tmp_CG_t3_",rand_string])]

    SeisCopy(d,r)
    push!(cost,InnerProduct(r,r))
    fp = open(cost_file,"w")
    write(fp,join(["began execution at: ",strftime(time()),"\n"]))
    write(fp,join([string(cost[1]),"\n"]))
    close(fp)
    LinearOperator(g,r,operators,parameters,adj=true)
    SeisCopy(g,s)
    SeisCopy(s,tt)
    CGStep(tt,s,a=[0.;0.;0.],b=[0.;0.;0.])
    SeisCopy(tt,m)
    
    gamma_old = InnerProduct(g,g)
    for iter = 1 : Niter    
        LinearOperator(s,t,operators,parameters,adj=false)
        CGStep(tt,s,a=[0.;0.;0.],b=mu)
        delta = InnerProduct(t,t) + InnerProduct(tt,tt)
        alpha = gamma_old/(delta + 1.e-20)
        CGStep(m,s,a=[1.;1.;1.],b=[alpha;alpha;alpha])  
        CGStep(r,t,a=[1.;1.;1.],b=-[alpha;alpha;alpha])
        CGStep(rr,tt,a=[1.;1.;1.],b=-[alpha;alpha;alpha])
        push!(cost,InnerProduct(r,r) + InnerProduct(rr,rr))
        fp = open(cost_file,"a")
        write(fp,join([string(cost[iter+1]),"\n"]))
        close(fp)
        LinearOperator(g,r,operators,parameters,adj=true)
        CGStep(g,rr,a=[1.;1.;1.],b=mu)
        gamma = InnerProduct(g,g)
        println("gamma=",gamma)
        beta = gamma/(gamma_old + 1.e-20)
        gamma_old = copy(gamma)
        CGStep(s,g,a=[beta;beta;beta],b=[1.;1.;1.])
    end
    SeisRemove(g);
    SeisRemove(s);
    SeisRemove(rr);
    SeisRemove(tt)
    SeisRemove(r);
    SeisRemove(t);

end


