using Seismic

pathin = "/Users/wenlei/Desktop/"   #change the directory to the veloctiy model
pathout= "/Users/wenlei/Desktop/"   #change the output directory

m = 1101 # number of rows of velocity model
n = 301  # number of column of velocity model
dt = 0.004 # sampling interval 0.004s
dx = 1     # grid size of velocity model 1m
fpp= 30    # doimant frequency

tmp_path = join([pathin "vp.bin"]) #path to your P-wave velocity file
fid = open(tmp_path, "r")
vp = read(fid, Float32, m*n)
vp = reshape(vp, m, n)
close(fid)

(dpp, rpp) = Synthetic(vp, dt, dx, fpp)
mP = [size(dpp,1), size(dpp,2)]
mR = [size(rpp,1), size(rpp,2)]

w  = Ricker(f0=fpp, dt=dt)                 # generate ricker wavelet
Wpp= convmtx(w, mR[1])                     # form convolution matrix
Wpp= sparse(kron(speye(mP[2]), Wpp))       # form convolution matrix for all traces


λ = 0.01  # the weight for sparse regularization term
β = 1     # the weight for horizontal continous term
niter = 50 # iteration numbers of fista

Dx = pardx(mR)  # second order derivative matrix along horiontal direction
Dx = sqrt(β)*Dx
H = vcat(Wpp, Dx)
b = vcat(vec(dpp), zeros(mR[1]*mR[2]))
param =["Hop"=>H, "nr"=>mR[1]*mR[2], "pow_niter"=>50]
println("Estimating maximum eigenvalue of operator by power method")
α = power_method(Hop, param)
println("sparse deconvolution solved by FISTA")
(rpp_decon, J) = fista(b, Hop, param, λ, α*0.98, niter)
rpp_decon = reshape(rpp_decon, mR[1], mR[2])

# plotting the result
vmax = maximum(dpp) * 0.6
vmin = -vmax
plotpar=["style"  =>"color" ,"wiggle_trace_increment" =>10 ,
         "cmap"   =>"seismic"  , "interpolation"=>"none"   ,
         "vmin"   => vmin   , "vmax" => vmax, "dpi" =>100  ,
         "aspect" =>"auto"  , "wbox" => 5   , "hbox"=>7  ,
         "title"  =>"", "name"=>"dpp.pdf",
         "xlabel" =>"trace" , "xunits"=>"number", "ox"=>1 , "dx"=>1,
         "ylabel" =>"Time", "yunits"=>"(s)", "oy"=>dt, "dy"=>dt,
         "titlesize"=>10, "labelsize"=>10, "ticksize"=>7     ]

# plotting the data
plotpar["name"] = join([pathout "dpp.pdf"])
SeisPlot(dpp, plotpar)

# plotting true PP-wave reflectivity
tmp = maximum(abs(rpp))
plotpar["vmax"] = tmp
plotpar["vmin"] = -tmp
plotpar["name"] = join([pathout "rpp_true.pdf"])
SeisPlot(rpp, plotpar)

# plotting the result of deconvolution
plotpar["title"] = "rpp_decon"
plotpar["titlesize"] = 20
plotpar["labelsize"] = 15
plotpar["ticksize"] = 15
plotpar["name"] = join([pathout "result_deconvolution.pdf"])
SeisPlot(rpp_decon, plotpar)
