# SeisWavelets.jl 
#
# A collection of functions to compute source wavelets and to estimate
# minimum phase  wavelets  via  Kolmogorov factorization
#                     
# Copyright (c) 2015 Signal Analysis and Imaging Group - University of Alberta
# http://saig.physics.ualberta.ca
#


module SeisWavelets 

function Ricker(param=Dict())
# 
# Ricker: Ricker wavelet of central frequency fc sampled every dt seconds
#
# IN   param={"dt"=>dt_value,"fc"=>fc_value}
# 
#      dt_value: sampling interval in secs
#      fc_value: central frequency [Hz]
# 
# OUT  t: time axis
#      w: the Ricker wavelet 
#
# Example with default values 
# 
#      (t,w) = SeisWavelets.Ricker()
#      plot(t,w) 
#      xlabel("Time (s)")
#
# Example with parameters provided by the user 
# 
#      param={"dt"=>0.004,"fc"=>20}
#      (t,w) = SeisWavelets.Ricker(param)
#      plot(t,w) 
#      xlabel("Time (s)")


	dt = get(param, "dt",0.002)
	fc = get(param, "fc",20)

	nw = 6./fc/dt
	nc = floor(nw/2)
	t = dt*[-nc:1:nc]
	b = (fc*pi*t).^2
	w = (1.-2.*b).*exp(-b)

 return t,w

end 

function Ormsby(param=Dict())
#
# Ormsby: Ormsby wavelet sampled every dt seconds with corner frequencies 
#         defined by the vector f=[f1,f2,f3,f4]. The final wavelet is multiply
#         by a Hamming window. 
#
# IN   param={"dt"=>dt_value,"f"=>f_value}
#
#      dt_value: sampling interval in secs
#      fc_value: a vector containing frequencies f1,f2,f3,f4  
#                defining the following amplitude response 
#
#      ^
#    1 |     ***************
#      |    *               *
#      |   *                 *
#      |  *                   * 
#      | *                     * 
#      -----------------------------> f 
#        f1  f2             f3  f4 
# 
# 
# OUT  t: time axis
#      w: the Ormsby wavelet 
#
# 
# Example with default parameters 
# 
#      (t,w) = Ormsby()
#      plot(t,w) 
#      xlabel("Time (s)")
#
# Example with parameters provided by the user 
# 
#      param={"dt"=>0.004,"f"=>[10.,20.,40.,50.]}
#      (t,w) = SeisWavelets.Ormsby(param)
#      plot(t,w) 
#      xlabel("Time (s)")

	dt = get(param,"dt",0.002)
	f = get(param,"f",[2.,10.,40.,60.])

	f1 = f[1]
	f2 = f[2]
	f3 = f[3]
	f4 = f[4]

	fc = (f2+f3)/2.
	nw = 6./(fc*dt)
	nc = floor(nw/2)
	t = dt*[-nc:1:nc]
        nw = 2*nc+1 
	a4 = (pi*f4)^2/(pi*(f4-f3))
	a3 = (pi*f3)^2/(pi*(f4-f3))
	a2 = (pi*f2)^2/(pi*(f2-f1))
	a1 = (pi*f1)^2/(pi*(f2-f1))

	u =   a4*(sinc(f4*t)).^2 - a3*(sinc(f3*t)).^2 
	v =   a2*(sinc(f2*t)).^2 - a1*(sinc(f1*t)).^2 

	w = u - v

        n = [0:1:nw-1]
        h = 0.54-0.46*cos(2*pi*n/(nw-1))
        wmax = maximum(w)
        w = w.*h/wmax

 return t,w

end

function Berlage(param=Dict())
#
# Berlage: Berlage wavelet
#
# IN   param={"dt"=>dt_value,"alpha"=>alpha_value,"m"=>m_value,"f0"=>f0_value,"phi0"=>phi0_value}
#
#      dt_value:    sampling interval in secs
#      alpha_value: a paramter of the Berlage wavelet
#      m_value:     another  parameter of the Berlage wavlet
#      f0_value:    central frequency [Hz]
#      phi0_value:  phase rotation [degrees]
# 
# OUT  t:  time axis
#      w:  the Berlage wavelet 
#
# Reference: David F. Aldridge, 1990, The Berlage wavelet: GEOPHYSICS, 55(11), 1508-1511.
#            doi: 10.1190/1.1442799
# 

	dt = get(param,"dt",0.002)
	alpha = get(param,"alpha",0.4)
	m = get(param,"m",2.)
	f0 = get(param,"f0",20.)
	phi0 = get(param,"phi0",3.14*90/180.)

	nw = 6./(f0*dt)
	t = dt*[0:1:nw-1]

	w = (t.^m).*exp(-alpha*t).*cos(2*pi*f0*t + 3.14*phi0/180.);

	wmax = maximum(w)
	w = w/wmax

 return t,w

end

function Kolmogoroff(param=Dict())
#
# kolmog: Kolmogoroff factorization. This program is used to transform a wavelet into its
#         minimum phase equivalent
#
# IN   param={"w"=>w_value}
#      w_value: a source wavelet with arbitraty phase given (a vector)
# 
# OUT  wmin:   the minimum phase wavelet  with amplitude spectrum identical to the 
#              amplitude spectrum of the input wavelet w 
#
# 
# Example with default values 
# 
#      (t,w) = SeismicWavelets.Ricker()
#       wmin = SeisimicWavelets.Kolmogoroff({"w"=>w})
#      plot(t,w,t,wmin)
#
# Reference: Jon F. Claerbout, 1976, Fundamentals of Geophysical Data Processing, McGraw Hill.

	w = get(param,"w",[])
	nw = length(w)
	nf = 8*nextpow2(nw)
	W = fft([w,zeros(nf-nw,1)])
	A = log(abs(W)+0.0000001)
	a = 2.*ifft(A)
	a[nf/2+2:nf] = 0.
	a[1]=a[1]/2.
	A = exp(fft(a))
	a = ifft(A)
	wmin = real(a[1:nw])

 return wmin

end
end
