function power_method(Hop, param)
  nr = get(param, "nr", 0)
  pow_niter = get(param, "pow_niter", 10)
  if nr == 0
    error("define the number of unknowns")
  end
  x = rand(nr)
  value = 0.
  n = 1.
  for k=1:pow_niter
    aux = Hop(x, param, 1)
    y = Hop(aux, param,-1)
    n = norm(x)
    x = y/n
    value = n
    println("iter $k, α=$value")
  end
  value = n
  return value
end
