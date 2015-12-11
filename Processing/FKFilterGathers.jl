function FKFilter(d,h,param)

nt = size(d,1)
nx = size(d,2)
dt = get(param,"dt",0.004)
dx = get(param,"dx",10)

pa = get(param,"pa",-.8*dx/dt)
pb = get(param,"pb",-.5*dx/dt)
pc = get(param,"pc",.5*dx/dt)
pd = get(param,"pd",.8*dx/dt)




m = fft(d,1)/sqrt(size(d,1))
m = fft(m,2)/sqrt(size(d,2))

nf = size(m,1)
dw = 2.*pi/nf/dt
nw = int(nf/2) + 1
fmax = get(param,"fmax",int(floor(0.5/dt)))
if(fmax*dt*nf < nw)
iw_max = int(floor(fmax*dt*nf))
else
iw_max = int(floor(0.5/dt))
end

nk = size(m,2)
dk = 2.*pi/nk/dx

for iw=2:iw_max+1
w = dw*(iw-1)
for ik=1:nk
ik < floor(nk/2) ? k = dk*ik : k = -(dk*nk - dk*ik)
p = k/w;
if (p<pa)
m[iw-1,ik] *= 0.
elseif (p >= pa && p < pb)
m[iw-1,ik] *= (p-pa)/(pb-pa)
elseif (p >= pb && p <= pc)
m[iw-1,ik] *= 1.
elseif (p > pc && p <= pd)
m[iw-1,ik] *= 1. - (p-pc)/(pd-pc)
else
m[iw-1,ik] *= 0.
end
end
end
fk_spec=abs(fftshift(m[1:fmax,:],2))
# symmetries
m = bfft(m,2)/sqrt(size(m,2))
for iw=nw+1:nf
m[iw,:] = conj(m[nf-iw+2,end:-1:1])
end
d = real(bfft(m,1)/sqrt(size(m,1)))
return d[1:nt,1:nx],fk_spec,h;
end

function FKFilterGathers(in,out,param)
SeisProcess(in,out,[FKFilter],param,"gather",["imx","imy"])
end