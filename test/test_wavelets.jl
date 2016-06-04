using PyPlot, Seismic

dt = 0.002
fn = 1/(2*dt)

# Ricker
w1 = Ricker(dt=dt)
nw = length(w1)
nc = floor(Int, nw/2)
t1 = dt*collect(-nc:1:nc)
nf = 8*nextpow2(nw)
df = 1/(nf*dt)
f1 = df*collect(0:1:nf-1)
wpad = cat(1,w1,zeros(nf-nw))
W1 = abs(fft(wpad))
W1 = W1/maximum(W1)

# Ormsby
w2 = Ormsby(dt=dt, f=[5.0,10.0,30.0,55.0])
nw = length(w2)
t2 = dt*collect(0:1:nw-1)
nf = 8*nextpow2(nw)
df = 1/(nf*dt)
f2 = df*collect(0:1:nf-1)
wpad = cat(1,w2,zeros(nf-nw))
W2 = abs(fft(wpad))
W2 = W2/maximum(W2)

# Berlage
w3 = Berlage(dt=dt)
nw = length(w3)
t3 = dt*collect(0:1:nw-1)
nf = 8*nextpow2(nw)
df = 1/(nf*dt)
f3 = df*collect(0:1:nf-1)
wpad = cat(1,w3,zeros(nf-nw))
W3 = abs(fft(wpad))
W3 = W3/maximum(W3)

figure("Suite of wavelets in time")

subplot(221)
plot(t1, w1)
axis("tight")
xlabel("Time (s)")
title("Ricker wavelet")

subplot(222)
plot(t2, w2)
axis("tight")
xlabel("Time (s)")
title("Ormsby wavelet")

subplot(223)
plot(t3, w3)
axis("tight")
xlabel("Time (s)")
title("Berlage wavelet")

tight_layout()

figure("Suite of amplitude spectrums")

subplot(221)
plot(f1, W1)
xlim([0, fn])
xlabel("Frequency (Hz)")
title("Ricker spectrum")
     
subplot(222)
plot(f2, W2)
xlim([0, fn])
xlabel("Frequency (Hz)")
title("Ormsby spectrum")

subplot(223)
plot(f3, W3)
xlim([0, fn])
xlabel("Frequency (Hz)")
title("Berlage spectrum")

#tight_layout()
