function Hop(in, param, iflag)

  H = get(param, "Hop", 0)
  if H==0
    error("define modeling operator")
  end

  if iflag == 1
    out = H * in
  else
    out = H'* in
  end
  return out
end
