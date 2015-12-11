function pardx(m)
#   m is the size of P-wave image
  Dx = spdiagm((ones(m[2]-1,1), -2*ones(m[2], 1), ones(m[2]-1,1)), [-1:1], m[2], m[2])
  Dx[[1, end]] = -Dx[[2, end-1]]
  Dx = sparse(kron(Dx, speye(m[1])))
  return Dx
end
