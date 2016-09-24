function SeisPSTM(m::Array{Float32,2},d::Array{Float32,2},param::Dict{Any,Any})

	v = get(param,"v",fill(1500.0f0,size(m)))
	adj = get(param,"adj",true)
	nsinc = int32(get(param,"nsinc",5))
	nsinc = isodd(nsinc) ? int32(nsinc) : int32(nsinc - 1) 
	aperture = float32(get(param,"aperture",1000.0))
	sx = get(param,"sx",fill(0.0f0,size(d,2)))
	sy = get(param,"sy",fill(0.0f0,size(d,2)))
	gx = get(param,"gx",fill(0.0f0,size(d,2)))
	gy = get(param,"gy",fill(0.0f0,size(d,2)))
	nt = int32(get(param,"nt",size(m,1)))
	ot = float32(get(param,"ot",0.0))
	dt = float32(get(param,"dt",0.001))
	nx = int32(get(param,"nx",size(m,2)))
	ox = float32(get(param,"ox",0.0))
	dx = float32(get(param,"dx",10.0))
	ny = int32(get(param,"ny",1))
	oy = float32(get(param,"oy",0.0))
	dy = float32(get(param,"dy",10.0))

	for itrace = 1 : size(d,2)
		if (adj == true)
			trace = d[:,itrace]
		end
		pstm_op(m,trace,v,sx[itrace],sy[itrace],gx[itrace],gy[itrace],nx,ox,dx,ny,oy,dy,nt,ot,dt,aperture,nsinc,adj)
		if (adj == false)
			d[:,itrace] = trace
		end
	end

end

function pstm_op(m::Array{Float32,2},trace::Array{Float32,1},v::Array{Float32,2},sx::Float32,sy::Float32,gx::Float32,gy::Float32,nx::Int32,ox::Float32,dx::Float32,ny::Int32,oy::Float32,dy::Float32,nt::Int32,ot::Float32,dt::Float32,aperture::Float32,nsinc::Int32,adj::Bool)

	trace = PhaseShift(trace,pi/2)
	trace1 = trace.*0
	trace2 = trace.*0
	Seismic.integrate(length(trace),trace,trace1,false)
	Seismic.integrate(length(trace),trace2,trace1,true)
	cmpx = (sx + gx)/2.0f0
	cmpy = (sy + gy)/2.0f0

	ixi = max(int32(1),int32((cmpx-aperture-ox)/dx))
	ixf = min(int32((cmpx+aperture-ox)/dx),nx)
	iyi = max(int32(1),int32((cmpy-aperture-oy)/dy))
	iyf = min(int32((cmpy+aperture-oy)/dy),ny)

	if (adj)
		for ix = ixi : ixf
			x = ox + ix*dx
			distx = abs(x - cmpx)
			for iy = iyi : iyf
				y = oy + iy*dy
				disty = abs(y - cmpy)
				dist = sqrt(distx^2 + disty^2)
				dists = sqrt((sx-x)^2 + (sy-y)^2)
				distg = sqrt((gx-x)^2 + (gy-y)^2)
				dists2 = dists*dists
				distg2 = distg*distg		
				for it = 1 : nt
					t0 = 0.5f0*(ot + it*dt)
					t02 = t0*t0
					v2 = v[it,(ix-1)*ny + iy]*v[it,(ix-1)*ny + iy]
					ts = sqrt(t02 + dists2/v2)
					tg = sqrt(t02 + distg2/v2)
					t = ts + tg
					geoms = sqrt(1.0f0/(t*v[it,(ix-1)*ny + iy]));
					obliq = sqrt(0.5f0*(1.0f0 + (t0*t0/(4.0f0*ts*tg)) - (1.0f0/(ts*tg))*sqrt(ts*ts - t0*t0/4.0f0)*sqrt(tg*tg - t0*t0/4.0f0)));
					fmax = 0.5f0*v2/dx/(dists/ts + distg/tg);
					jt_float = (t-ot)/dt
					jt = int(jt_float)
					k = max(round(1/fmax/dt-1),0)
					it1 = jt-int(k)-1
					it2 = jt+int(k)+1
					if nsinc < it1 && it2 < nt - nsinc
						for isinc = - (nsinc - 1)/2 : (nsinc - 1)/2
							m[it,(ix-1)*ny + iy] += -geoms*obliq*sinc(float32( (isinc)/nsinc ))*trace2[it1 + isinc]/((k+1)^2)
							m[it,(ix-1)*ny + iy] += 2*geoms*obliq*sinc(float32( (isinc)/nsinc ))*trace2[jt + isinc]/((k+1)^2)
							m[it,(ix-1)*ny + iy] += -geoms*obliq*sinc(float32( (isinc)/nsinc ))*trace2[it2 + isinc]/((k+1)^2)
						end
						#m[it,(ix-1)*ny + iy] += 2*geoms*obliq*trace[jt]
					end
				end
			end
		end
	else
		error("write it now")
	end

end
