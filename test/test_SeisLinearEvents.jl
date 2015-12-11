using PyPlot, Seismic

# generate synthetic 2d data consisting of linear events
# ======================================================
param = ["nt"=>500,"nx1"=>500,
	 "tau1"=>[0.4 1.0],"tau2"=>[0. 0.],"tau3"=>[0. 0.],"tau4"=>[0. 0.],
	 "v1"=>[3500. -4000],"v2"=>[99999. 99999.],"v3"=>[99999. 99999.],"v4"=>[99999. 99999.],
     "amp"=>[1. -0.5], "f0"=>[20. 20.]];
d = SeisLinearEvents(param);

plotpar = {"style"=>"overlay",
           "wiggle_trace_increment"=>10,
           "xcur"=>1.4,
           "vmin"=>-2,"vmax"=>2,
           "aspect"=>"auto",
           "xlabel"=>"X","xunits"=>"meters","ox"=>0,"dx"=>10,
           "ylabel"=>"Time","yunits"=>"seconds","oy"=>0,"dy"=>0.004,
           "wbox"=>8,"hbox"=>8,
           "cmap"=>"seismic"};

d_2d = squeeze(d,[5 4 3]);

plotpar["style"]="color"; plotpar["title"]="color";  
SeisPlot(d_2d,plotpar);
plotpar["style"]="wiggles"; plotpar["title"]="wiggles";  
SeisPlot(d_2d,plotpar);
plotpar["style"]="overlay"; plotpar["title"]="overlay";  
SeisPlot(d_2d,plotpar);

