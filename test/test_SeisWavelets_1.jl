
# Compute and display seismic wavelets 

using PyPlot,Seismic

      figure("Suite of wavelets")

      param={"dt"=>0.004,"fc"=>20}
      (t,w) = SeisWavelets.Ricker(param)

      subplot(221)
      plot(t,w)
      axis("tight")
      xlabel("Time (s)")
      title("Ricker")

      param={"dt"=>0.002,"f"=>[5.,10.,30.,55.]}
      (t,w) = SeisWavelets.Ormsby(param)

      subplot(222)
      plot(t,w)
      axis("tight")
      xlabel("Time (s)")
      title("Ormsby")

      param={"dt"=>0.004,"f0"=>20, "phi0"=>90., "alpha"=>50., "m"=>2.}
      (t,w) = SeisWavelets.Berlage(param)

      subplot(223)
      plot(t,w)
      axis("tight")
      xlabel("Time (s)")
      title("Berlage")

