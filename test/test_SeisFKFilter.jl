using Seismic, PyCall
@pyimport matplotlib.pyplot as pl

#Generate seismic data

nt = 256; 
nx = convert(Int64,128);
h = Array{Header}(nx);
for ix = 1:nx
  h[ix] = Seismic.InitSeisHeader();
  h[ix].tracenum = ix;
  h[ix].d1 = 0.004;
  h[ix].n1 = nt;
  h[ix].imx = ix;
  h[ix].imy = 0;
  h[ix].ihx = 0;
  h[ix].ihy = 0;
end

dt=0.004;
dx=5;
f_niquest=Int(floor(1/dt))

Tmax=(nt-1)*dt

data,ext=SeisLinearEvents(nt=nt,dt=dt,nx1=nx,dx1=dx,tau=[Tmax/4., Tmax/3])#,p1=[0.0001])

#Transform to FK domain
m=fft(data,1)/sqrt(size(data,1));
m=fft(m,2)/sqrt(size(data,2))
m=fftshift(m)

m=m[Int(floor(f_niquest/2)):f_niquest,:]

# apply FK filter
data_fk= SeisFKFilter(data;dt=dt,dx=dx,va=-1500,vb=-2500,vc=2500,vd=1500)

#Transform back to TX domain.

m=fft(data_fk,1)/sqrt(size(data_fk,1));
m=fft(m,2)/sqrt(size(data_fk,2))
m=fftshift(m)

m=m[Int(floor(f_niquest/2)):f_niquest,:]

#>>>>>>>>>>>>>>>Plot<<<<<<<<<<<<<<<<

f1 = pl.figure(1)
pl.imshow(data[1:nt,:],cmap="PuOr",aspect="auto",vmin=-1,vmax=1)
pl.title("data before fk filter")
pl.xlabel("X (m)")
pl.ylabel("T (s)")
pl.show()



f2 = pl.figure(2)
pl.imshow(abs.(m),cmap="PuOr",aspect="auto",vmin=-1,vmax=1)
pl.title("fk spectrum before filter")
pl.xlabel("K (1/m)")
pl.ylabel("W (1/s)")
pl.show()


f3 = pl.figure(3)
pl.imshow(data_fk[1:nt,:],cmap="PuOr",aspect="auto",vmin=-1,vmax=1)
pl.title("data after fk filter")
pl.xlabel("X (m)")
pl.ylabel("T (s)")
pl.show()



f4 = pl.figure(4)
pl.imshow(abs.(m),cmap="PuOr",aspect="auto",vmin=-1,vmax=1)
pl.title("fk spectrum after filter")
pl.xlabel("K (1/m)")
pl.ylabel("W (1/s)")
pl.show()

pl.show()
