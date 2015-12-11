
# Compute and display  a ricker wavelet and compute its minimum phase
# equivalent via Kolmogoroff factorization

using PyPlot,Seismic

      figure("Mimimum Phase test")

      param={"dt"=>0.004,"fc"=>20}
      (t,w) = SeisWavelets.Ricker(param)

      subplot(221)
      plot(t,w)
      axis("tight")
      xlabel("Time (s)")
      title("Ricker")

      param={"w"=>w}
      wmin = SeisWavelets.Kolmogoroff(param)

      subplot(222)
      t = [0:1:length(wmin)-1]*0.004
      plot(t,wmin)
      axis("tight")
      xlabel("Time (s)")
      title("Minimum Phase equivalent")
