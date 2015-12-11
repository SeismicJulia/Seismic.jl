
using PyPlot, Seismic

# Test multiple supression via the parabolic Radon Transform 
# 
# 1-  Make a Radon Gather 
# 2-  Synthetize an NMO-corrected gather using the forward Radon transform 
# 3-  Use the NMO-corrected gather to recover the Radon panel via LS inversion
# 4-  Mute primeries in the Radon gather
# 5-  Obtain multiples by forward modeliing 4
# 6-  Substrat multiples from data to obtain primaries 
#
# References:
#
#  Hampson, D., 1986, Inverse velocity stacking for multiple elimination: Canadian Journal of Exploration Geophysics, 22, 44-55.
#  Sacchi, M.D. and Ulrych, T.J., 1995, High-resolution velocity gathers and offset space reconstruction: Geophysics, 60, 1169-1177.
#

#### To do -- add HR version of PRT to program SeisRadonInverse ######


        close("all")

# 1- Make an ideal Radon gather 

        t,source=SeisWavelets.Ricker()
        ns = length(source)
        np = 100
        nt = 800
        nw = nt
        nh = 80
        dt = 4./1000
        h = linspace(0.,1000.0,nh)
        p = linspace(-0.04,2.2,np)
        m = zeros(Float64,nt,np)
        for is = 1:ns
         m[160+is,3 ]=source[is]
         m[ 60+is,44]=-source[is]
         m[280+is,44]=-source[is]
         m[240+is, 3]=-source[is]
        end
        href = maximum(abs(h))

# 2- Model the data 

# Pass arguments of Radon forward transform via dictionary 

    param = {"dt"=>dt,"p"=>p, "h"=>h, "href"=>href, "fhigh"=>40., "flow"=>3., "parab"=>true} 

    d = SeisRadonForward(m,param);

# Add trade-off parameter for the inversion to dictionary 

    param["mu"]=0.001

# 3-  Recover Radon gather via inversion

    m = SeisRadonInverse(d,param);

# 4- Filter primaries and keep multiples in Radon gather  

    mf = copy(m);
    mf[:,1:16]=0.

# 5- Model multiples 

    d_mult = SeisRadonForward(mf,param);
 
# 6- Compute primaries

    d_prim = d-d_mult 

    dp = p[2]-p[1]
    dh = h[2]-h[1]

# Plotting  

    Plot_param = {"title"=>"Data", "xlabel"=>"Offset (m)", "ylabel"=>"Time (s)", "ox"=>h[1], "dx"=>dh, "dy"=>dt}
    SeisPlot(d, Plot_param)

    Plot_param["title"]="Multiples"
    SeisPlot(d_mult, Plot_param)

    Plot_param["title"]="Primaries"
    SeisPlot(d_prim, Plot_param)

    Plot_param = {"title"=>"Radon Gather", "xlabel"=>"Curvature (s)", "ylabel"=>"Time (s)", "ox"=>p[1], "dx"=>dp, "dy"=>dt}
    SeisPlot(m, Plot_param)
