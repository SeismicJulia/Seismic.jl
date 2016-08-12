function ConjugateGradients(d,operators,parameters;Niter=10,mu=0,tol=1.0e-15)
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	cost = Float64[]
	r = copy(d)
	g = LinearOperator(r,operators,parameters,adj=true)
	m = zeros(g)
	s = copy(g)
	gamma = InnerProduct(g,g)
	gamma00 = gamma
    cost0 = InnerProduct(r,r)
	push!(cost,1.0)
	for iter = 1 : Niter
		t = LinearOperator(s,operators,parameters,adj=false)
		delta = InnerProduct(t,t) + mu*InnerProduct(s,s)		
		if delta <= tol 
			println("delta reached tolerance, ending at iteration ",iter)
			break;
		end		
		alpha = gamma/delta
		m = m + alpha*s
		r = r - alpha*t
		g = LinearOperator(r,operators,parameters,adj=true)
		g = g - mu*m
		gamma0 = copy(gamma)
		gamma = InnerProduct(g,g)
        cost1 = InnerProduct(r,r) + mu*InnerProduct(m,m)
        push!(cost,cost1/cost0)
		beta = gamma/gamma0
		s = beta*s + g
		if (sqrt(gamma) <= sqrt(gamma00) * tol)
			println("tolerance reached, ending at iteration ",iter)
			break;
		end
	end

	return m, cost
end

function ConjugateGradients(m::ASCIIString,d::ASCIIString,operators,parameters,cost_file::ASCIIString;Niter=10,mu=0,tol=1.0e-15)
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	cost = Float64[]
	rand_string = string(round(Int,rand()*100000))
	g = join(["tmp_CG_g_",rand_string])
	s = join(["tmp_CG_s_",rand_string])
	r = join(["tmp_CG_r_",rand_string])
	t = join(["tmp_CG_t_",rand_string])
	SeisCopy(d,r)
	fp = open(cost_file,"w")
	@compat write(fp,join(["began execution at: ",Libc.strftime(time()),"\n"]))
	close(fp)
	LinearOperator(g,r,operators,parameters,adj=true)
	SeisCopy(g,s)
	SeisCopy(g,m)
	CGStep(m,g,a=0.0,b=0.0)
	gamma = InnerProduct(g,g)
	gamma00 = gamma
    cost0 = InnerProduct(r,r)
	push!(cost,1.0)
	fp = open(cost_file,"a")
	write(fp,join([string(cost[1]),"\n"]))
	close(fp)
	for iter = 1 : Niter	
		LinearOperator(s,t,operators,parameters,adj=false)
		delta = InnerProduct(t,t) + mu*InnerProduct(s,s)		
		alpha = gamma/delta
		CGStep(m,s,a=1.0,b=alpha)  
		CGStep(r,t,a=1.0,b=-alpha)
		LinearOperator(g,r,operators,parameters,adj=true)
		CGStep(g,m,a=1.0,b=-mu)
		gamma0 = copy(gamma)
		gamma = InnerProduct(g,g)
        cost1 = InnerProduct(r,r) + mu*InnerProduct(m,m)
        push!(cost,cost1/cost0)
		fp = open(cost_file,"a")
		write(fp,join([string(cost[iter+1]),"\n"]))
		close(fp)
		beta = gamma/gamma0
		CGStep(s,g,a=beta,b=1.0)
		if (sqrt(gamma) <= sqrt(gamma00) * tol)
			println("tolerance reached, ending at iteration ",iter)
			break;
		end
	end
	SeisRemove(g)
	SeisRemove(s)
	SeisRemove(r)
	SeisRemove(t)
	
end

function ConjugateGradients(m::Array{ASCIIString,1},d::Array{ASCIIString,1},operators,parameters,cost_file::ASCIIString;Niter=10,mu=[0.0;0.0],tol=1.0e-15)
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	cost = Float64[]
	rand_string = string(round(Int,rand()*100000))
	g = [join(["tmp_CG_g1_",rand_string]);join(["tmp_CG_g2_",rand_string])]
	s = [join(["tmp_CG_s1_",rand_string]);join(["tmp_CG_s2_",rand_string])]
	r = [join(["tmp_CG_r1_",rand_string]);join(["tmp_CG_r2_",rand_string])]
	t = [join(["tmp_CG_t1_",rand_string]);join(["tmp_CG_t2_",rand_string])]
	SeisCopy(d,r)
	fp = open(cost_file,"w")
	@compat write(fp,join(["began execution at: ",Libc.strftime(time()),"\n"]))
	close(fp)
	LinearOperator(g,r,operators,parameters,adj=true)
	SeisCopy(g,s)
	SeisCopy(g,m)
	CGStep(m,g,a=[0.0;0.0],b=[0.0;0.0])
	gamma = InnerProduct(g,g)
	gamma00 = gamma
    cost0 = InnerProduct(r,r)
	push!(cost,1.0)
	fp = open(cost_file,"a")
	write(fp,join([string(cost[1]),"\n"]))
	close(fp)
	for iter = 1 : Niter	
		LinearOperator(s,t,operators,parameters,adj=false)
		delta = InnerProduct(t,t) + mu[1]*InnerProduct(s[1],s[1]) + mu[2]*InnerProduct(s[2],s[2])	
		alpha = gamma/delta
		CGStep(m,s,a=[1.0;1.0],b=[alpha;alpha]) 
		CGStep(r,t,a=[1.0;1.0],b=[-alpha;-alpha])
		LinearOperator(g,r,operators,parameters,adj=true)
		CGStep(g,m,a=[1.0;1.0],b=[-mu[1];-mu[2]])
		gamma0 = copy(gamma)
		gamma = InnerProduct(g,g)
        cost1 = InnerProduct(r,r) + mu[1]*InnerProduct(m[1],m[1]) + mu[2]*InnerProduct(m[2],m[2])
        push!(cost,cost1/cost0)
		fp = open(cost_file,"a")
		write(fp,join([string(cost[iter+1]),"\n"]))
		close(fp)
		beta = gamma/gamma0
		CGStep(s,g,a=[beta;beta],b=[1.0;1.0])
		if (sqrt(gamma) <= sqrt(gamma00) * tol)
			println("tolerance reached, ending at iteration ",iter)
			break;
		end
	end
	SeisRemove(g);
	SeisRemove(s);
	SeisRemove(r);
	SeisRemove(t);
	
end
