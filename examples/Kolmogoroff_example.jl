using PyPlot, Seismic

dt = 0.002
w = Ricker(dt=dt)
wmin = SeisKolmogoroff(w)

nw = length(w)
t = dt*collect(0:1:nw-1)

figure("Kolmogoroff factorization", figsize=(6, 6))

subplot(211)
plot(t, w)
axis("tight")
xlabel("Time (s)")
title("Ricker wavelet")

subplot(212)
plot(t,wmin)
axis("tight")
xlabel("Time (s)")
title("Minimum Phase equivalent")

tight_layout()
