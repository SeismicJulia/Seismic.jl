"""
    lambda = power_method(m, operators, parameters, max_iter=50)

random algorithm to estimate the max eigenvalue of linear operators

# Arguments
* `m`: random model parameters
* `operators  :: Array{Function,1}`: array of linear operators
* `parameters :: Array{Dict,1}`: array of dictionary
* the length of parameters must be equal to that of operators, the ith element of
* parameters contains the argument to ith linear operator in operators.

# Keywords arguments
* `max_iter :: Int64` : number of iterations

# Output
* `lambda :: Float64`: max eigenvalue
"""
function power_method(m, operators, parameters; max_iter=50)
    lambda = 0.0
    for iter = 1 : max_iter
        d = LinearOperator(m, operators, parameters, adj=false)
        m1= LinearOperator(d, operators, parameters, adj=true)
        lambda = norm(m)
        m = m1 / lambda
        println("iteration: $iter, maximum eig: $lambda")
    end
    return lambda
end

function softThresh(m, t::Float64)
    tmp = abs(m) - t
    tmp = (tmp + abs(tmp)) / 2
    m   = sign(m) .* tmp
    return m
end

"""
    (m, J) = FISTA(d, operators, parameters, mu, lambda, max_iter=50)

FISTA solver for ||d - op1 * op2 * op3 * m||2^2 + mu||m||1

# Arguments
* `d`: measured data
* `operators  :: Array{Function,1}`: array of linear operators
* `parameters :: Array{Dict,1}`    : array of dictionary
* the length of parameters must be equal to that of operators, the ith element of
* parameters contains the argument to ith linear operator in operators.
* `mu :: Float64`: measured data
* `lambda :: Float64`: max eigenvalue estimated by power's method.

# Keywords arguments
* `max_iter :: Int64` : number of iterations

# Output
* `m `: estimated model parameters
* `J :: Array{Float64,1}`: the value of objective function at each iteration
"""
function FISTA(d, operators, parameters, mu::Float64, lambda::Float64; max_iter=50)
    J = zeros(max_iter+1)
    T = mu / (2*lambda)
    cost = dot(vec(d), vec(d)); J[1] = cost;
    println("iterations: 0, object value $cost")
    m = LinearOperator(-d, operators, parameters, adj=true)
    m0= zeros(m); yk = zeros(m)
    m = - m/lambda
    m = softThresh(m, T)
    t = 1.0; t0 = t;
    t = (1 + sqrt(1+4*t^2)) / 2
    yk= m
    for k = 1 : max_iter
        m0[:] = m[:]
        df = LinearOperator(m, operators, parameters, adj=false)
        df = df - d
        cost = dot(vec(df), vec(df)) + mu*sum(abs(m)); J[k+1] = cost;
        println("iterations: $k, object value $cost")
        m  = LinearOperator(df, operators, parameters, adj=true )
        m  = yk - m/lambda
        m  = softThresh(m, T)
        t0 = t
        t  = (1 + sqrt(1+4*t^2)) / 2
        yk = m + (t0-1)/t * (m-m0)
    end
    return m, J
end
