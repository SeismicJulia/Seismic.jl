function SeisMWNI(in,h::Array{Header,1},param::Dict{Any,Any})

	style = get(param,"style","sxsygxgy")
	padt = get(param,"padt",2)
	padx = get(param,"padx",2)
	param["operators"] = [ApplyDataWeights fft_op ApplyModelWeights ]
	nt = size(in,1)
	nx = size(in,2)
	ot = h[1].o1
	dt = h[1].d1
	nf = padt*nextpow2(nt)
	dw = 2.*pi/nf/dt
	nw = int(nf/2) + 1
	fmax = get(param,"fmax",int(floor(0.5/dt)))
	if(fmax*dt*nf < nw) 
		iw_max = int(floor(fmax*dt*nf))
	else 
		iw_max = int(floor(0.5/dt))
	end

	if (style == "sxsygxgy")
		key = ["t","isx","isy","igx","igy"]
	elseif (style=="mxmyhxhy")
		key = ["t","imx","imy","ihx","ihy"]
	elseif (style=="mxmyhaz")
		key = ["t","imx","imy","ih","iaz"]
	elseif (style=="sxsyhxhy")
		key = ["t","isx","isy","ihx","ihy"]
	elseif (style=="gxgyhxhy")
		key = ["t","igx","igy","ihx","ihy"]
	elseif (style=="sxsyhaz")
		key = ["t","isx","isy","ih","iaz"]
	elseif (style=="gxgyhaz")
		key = ["t","igx","igy","ih","iaz"]
	else
		error("style not defined.")
	end
	min_ix1 = getfield(h[1],symbol(key[2]))
	max_ix1 = getfield(h[nx],symbol(key[2]))
	nx1 = max_ix1 - min_ix1 + 1
	min_ix2 = getfield(h[1],symbol(key[3]))
	max_ix2 = getfield(h[nx],symbol(key[3]))
	nx2 = max_ix2 - min_ix2 + 1
	min_ix3 = getfield(h[1],symbol(key[4]))
	max_ix3 = getfield(h[nx],symbol(key[4]))
	nx3 = max_ix3 - min_ix3 + 1
	min_ix4 = getfield(h[1],symbol(key[5]))
	max_ix4 = getfield(h[nx],symbol(key[5]))
	nx4 = max_ix4 - min_ix4 + 1
	d = reshape(in,nt,nx1,nx2,nx3,nx4)
	nx1 > 1 ? nk1 = padx*nextpow2(nx1) : nk1 = 1
	nx2 > 1 ? nk2 = padx*nextpow2(nx2) : nk2 = 1
	nx3 > 1 ? nk3 = padx*nextpow2(nx3) : nk3 = 1
	nx4 > 1 ? nk4 = padx*nextpow2(nx4) : nk4 = 1
	nk = nk1*nk2*nk3*nk4
	# generate sampling operator from the padded data
	T,h_tmp = CalculateSampling(d)
	T = T[1,:,:,:,:]
	if (sum(T[:])/length(T[:]) < 0.05)
		println(sum(T[:])/length(T[:]))
		return in,h
	else
		T = Pad5D(T,1,nk1,nk2,nk3,nk4)
		T = squeeze(T[1,:,:,:,:],1)
		param["wd"] = T
		d = Pad5D(d,nf,nk1,nk2,nk3,nk4)
		D = fft(d,1)
		for iw=1:iw_max
			#println(iw,"/",iw_max)
			x = squeeze(D[iw,:,:,:,:],1)
			y = copy(x)
			param["adj"] = true
			Y = fft_op(y,param)
			#param["wm"] = 1.
			Y,misfit = IRLS(Y.*0,x,param)
			param["adj"] = false
			y = fft_op(Y,param)
			D[iw,:,:,:,:] = y
		end
		# symmetries
		for iw=nw+1:nf
			D[iw,:,:,:,:] = conj(D[nf-iw+2,:,:,:,:])
		end 
		d = ifft(D,1)
		d = real(d[1:nt,1:nx1,1:nx2,1:nx3,1:nx4])
		out = reshape(d,nt,nx1*nx2*nx3*nx4)
		return out,h
	end
end
