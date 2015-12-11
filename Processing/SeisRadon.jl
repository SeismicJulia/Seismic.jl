function SeisRadonForward(d,h,param=Dict()) 
	(nt,np) = size(m)
        M = fft(m,1)
        D = zeros(Complex64,nw,nh)



        for iw = iw_low:iw_high

         w  = 2.*pi*(iw-1)/(nw*dt)

         L = zeros(Complex64,nh,np)       # Forward operator 


          for ip = 1:np
          for ih = 1:nh
                phi = w*p[ip]*(h[ih]/href)^2
      	   	L[ih,ip] = exp(-im*phi)
          end
          end 
	
         x = M[iw,:].'
         y = L*x

         D[iw,:] = y.'

        end 

          for iw = int(nw/2)+2:nw
           D[iw,:] = conj(D[nw-iw+2,:])
          end

         d = real(ifft(D,1))

# Go Back ...

	(nt,nh) = size(d)
        D = fft(d,1)
        M = zeros(Complex64,nw,np)


        for iw = iw_low:iw_high

         w  = 2.*pi*(iw-1)/(nw*dt)

         L = zeros(Complex64,nh,np)       # Forward operator 


          for ip = 1:np
          for ih = 1:nh
                phi = w*p[ip]*(h[ih]/href)^2
      	   	L[ih,ip] = exp(-im*phi)
          end
          end 

         R = L'*L; 
         r0 = trace(R)/np
         mu = 0.01
         y = D[iw,:].'
         xa = L'*y
         Q = mu*r0*eye(np)
          for iter = 1:4 
         x  = (R + Q)\xa
         q = 1./ ( abs(x[:,1]).^2 + 0.001)
         Q = diagm(q)
         end
         M[iw,:] = x.'

        end 
          for iw = int(nw/2)+2:nw
           M[iw,:] = conj(M[nw-iw+2,:])
          end

         m = real(ifft(M,1))



 SeisPlot(d)
 SeisPlot(m)

return
 
return
 
