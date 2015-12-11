function SeisRadonInverse(d,param=Dict()) 
#
# SeisRadonInverse: frequency domain forward Radon Transform
#
# *** This programs transform a CMP Gather (t-h)  to tau-p using least-squares inversion *** 
#
#
# IN   d:        data,  d[1:nt,1:nh] nt is number of samples and nh is number of offsets 
#
#      param={"dt"=>dt,"parab"=>true/false, "p"=>p, "h"=>h, "flow"=>flow, "fhigh"=>fhigh, "href"=>href}
#
#      dt:       sampling interval in secs
#      parab:    true --> parabolic transform, false --> linear transform 
#      h:        offset vector, h[1:nh]
#      p:        if parab=true:  vector of residual moveout at reference offset in seconds (also called curvature)  
#                if parab=false: vector of ray parameters  reference (s/m)  
#                p[1:np] np is number of curvatures or ray parameters
#      flow:     min frequency in the data (Hz)
#      fhigh:    max frequency in the data (Hz)
#      href:     reference offset for parabolic Radon Transform 
#      mu:       trade off parameter or damping for the least squares inversio
#
# OUT  m:        the Radon panel of size m[1:nt, 1:np]
#

        parab = get(param,"parab",true)   # default is parabolic
        href = get(param,"href",1.)

        if param["parab"] == true
         I = 2                          # Parabolic moveout
        else
         I = 1                          # Linear moveout
         href = 1
        end


        dt = get(param,"dt",0.002)
        p = get(param,"p",[])
        h = get(param,"h",[])
        flow = get(param,"flow",2.)
        fhigh = get(param,"fhigh",40.)
        mu = get(param,"mu",0.001);

        (nt,nh) = size(d)
        nw = nextpow2(nt)
        np = length(p)
        d = [d;zeros(nw-nt,nh)]

        D = fft(d,1)

        iw_low = int(flow*dt*nw+1)
        iw_high = int(fhigh*dt*nw+1)
        M = zeros(Complex64,nw,np)

        for iw = iw_low:iw_high

         w  = 2.*pi*(iw-1)/(nw*dt)

         L = zeros(Complex64,nh,np)       # Forward operator 

          for ip = 1:np
          for ih = 1:nh
                phi = w*p[ip]*(h[ih]/href)^I
      	   	L[ih,ip] = exp(-im*phi)
          end
          end 

         R = L'*L; 
         r0 = trace(R)/np
         y = D[iw,:].'
         xa = L'*y
         Q = mu*r0*eye(np)
         x  = (R + Q)\xa
         M[iw,:] = x.'

        end 
          for iw = int(nw/2)+2:nw
           M[iw,:] = conj(M[nw-iw+2,:])
          end

         m = real(ifft(M,1))

         m = m[1:nt,:];

return m
end 
 
