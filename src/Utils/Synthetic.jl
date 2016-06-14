function Synthetic(vp, dt, dx, fpp)

#   the dimension of velocity model
  (m,n) = size(vp)
  nr = m-1
  rc_pp = zeros(nr,n)

#  compute reflection coefficients
  for in = 1:n
     for k = 1:nr
       rc_pp[k, in] = (vp[k+1,in]-vp[k,in]) / (vp[k+1,in]+vp[k,in])
     end
  end

  tpp = zeros(nr,n)

#   the time location of reflector
  for in = 1:n
    for k = 1:nr
      for l =1: k
        tpp[k, in] = tpp[k, in] + 2*dx/vp[l, in]
      end
    end
  end

  ntp = int(floor( maximum(tpp)/dt)) + 30
  rpp = zeros(ntp,n)

#   put reflection coefficient to right time location
  for in = 1:n
    for k = 1:nr
      jp = int(floor(tpp[k, in]/dt)+1)
      if abs(rpp[jp, in]) < abs(rc_pp[k, in])
         rpp[jp,in] = rc_pp[k, in]
      end
    end
  end

#   convolve rpp with wavelet
  wp = Ricker(f0=fpp, dt=dt)
  lenw = length(wp)
  lenrp= size(rpp, 1)
  len  = lenw+lenrp-1
  println("dpp length $len")
  dpp  = zeros(len, n)

#   convolve rpp with wavelet
  for in = 1:n
    dpp[:, in] = conv(rpp[:, in], wp)
  end
  return dpp, rpp
end
