"""
**SeisRadonForward**

*Transform a tau-p gather to a time-offset gather using a frequency domain forward Radon operator*

**IN**   

* m: Radon panel,  m[1:nt,1:nh] nt is number of samples and np is number of curvatures or ray parameters
* param={"dt"=>dt,"parab"=>true/false, "p"=>p, "h"=>h, "flow"=>flow, "fhigh"=>fhigh, "href"=>href}
  * dt: sampling interval in secs
  * parab: true --> parabolic transform, false --> linear transform
  * h: offset vector, h[1:nh]
  * if parab=true:  p is a vector of residual moveout at reference offset in seconds (also called curvature)
  * if parab=false: P is a vector of ray parameters  reference (s/m), p[1:np] where np is number of curvatures or ray parameters
  * flow: min frequency in the data (Hz)
  * fhigh: max frequency in the data (Hz)
  * href: reference offset for parabolic Radon Transform

**OUT**  

* d: data synthetized via forward Radon modeling, d[1:nt, 1:nh]
"""
function SeisRadonForward(m,param=Dict()) 

        parab = get(param,"parab",true)   # default is parabolic 
        href = get(param,"href",1000.) 

        if param["parab"] == true 
         I = 2                          # Parabolic moveout 
        else
         I = 1                          # Linear moveout 
         href = 1.
        end

        dt = get(param,"dt",0.002)
        p = get(param,"p",[])
        h = get(param,"h",[])
        flow = get(param,"flow",2.)
        fhigh = get(param,"fhigh",40.)

	(nt,np) = size(m)
        nw = nextpow2(nt)
        nh = length(h)
        m = [m;zeros(nw-nt,np)]

        M = fft(m,1)
        
        iw_low = int(flow*dt*nw+1)
        iw_high = int(fhigh*dt*nw+1)
        D = zeros(Complex64,nw,nh)

        for iw = iw_low:iw_high

         w  = 2.*pi*(iw-1)/(nw*dt)

         L = zeros(Complex64,nh,np);      # Forward operator 


          for ip = 1:np
          for ih = 1:nh
                phi = w*p[ip]*(h[ih]/href)^I
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
         d = d[1:nt,:];

   return d
 
end
