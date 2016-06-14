using PyPlot,Seismic
d,e = SeisLinearEvents();
SeisPlot(d, style="wiggles", xcur=1, wiggle_trace_increment=5, oy=0, dy=0.004);
