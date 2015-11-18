function fista(y, Hop, param, λ, α, niter)
#   minimize J(x) = || y-Hx||₂² + λ|x|₁
#   input:
#      y observed data
#      Hop operator and its adjoint
#      param all stuff needed by operator
#      α  α>= max(eig(H'*H))
#      niter number of iterations
#   output:
#        x - result of inversion
#        J - objective function
  J = zeros(niter)
  x = zeros(Hop(y, param, -1))
  T = λ / (2*α)

  t = 1
  yk= x
  for k=1:niter
    tmpx = x
    Hx   = Hop(yk, param, 1)
    x    = thresh(yk+(Hop(y-Hx, param, -1))/α, "s", T)
    J[k] = sum(abs(Hx[:]-y[:]).^2) + λ*sum(abs(x[:]))
    tmp  = J[k]
    println("iter $k, J=$tmp")

    tmpt = t
    t    = (1+sqrt(1+4*t^2))/2
    yk   = x + (tmpt-1)/t*(x-tmpx)
  end

  return x, J
end
