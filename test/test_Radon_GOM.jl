using PyPlot, Seismic

close("all")

# Multiple supresssion using the Parabolic Radon transform 

# Download data set from web site (1 NMO-corrected CMP gather from the GOM)

    download("http://seismic.physics.ualberta.ca/data/gom_cdp_nmo.su","gom_cdp_nmo.su")

# Convert SU to internal SeismicJulia format

    SegyToSeis("gom_cdp_nmo.su","gom_cdp_nmo",["format"=>"su","input_type"=>"ieee","swap_bytes"=>false])

# Read data d and trace headers h

    (d,h) = SeisRead("gom_cdp_nmo")

     d = d[750:1500,:]

# Extracct offset  from headers

   offset=Seismic.ExtractHeader(h,"h")
   offset = offset*0.3048                  # I like SI units 

   nt,nh = size(d)

   href = maximum(abs(offset))
   dt = h[1].d1

# Curvatures in seconds (residual moveout at reference offset)

   np = 200
   p = linspace(-0.2,0.8,np)

   param = {"dt"=>0.004,"p"=>p, "h"=>offset, "href"=>href, "fhigh"=>40., "flow"=>3., "parab"=>true} 
   param["mu"]=0.1

   m = SeisRadonInverse(d,param);

   plot_param = {"ox"=>p[1], "dx"=>p[2]-p[1], "dy"=>dt}
   SeisPlot(m,plot_param)

   plot_param = {"ox"=>offset[1], "dx"=>offset[2]-offset[1], "dy"=>dt}
   SeisPlot(d,plot_param)

# Apply mutes

   mf = copy(m)
   mf[:,1:70]=0

# Forward model multiples from Radon panel 

   d_mult = SeisRadonForward(mf,param);
   

# Compute and plot primaries 

   d_prim = d - d_mult


   SeisPlot(d_prim,plot_param)
