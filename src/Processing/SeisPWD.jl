function SeisPWD(in;w1=10,w2=10,dz=1,dx=1)
	# dip estimation by Plane Wave Destruction
	# see http://sepwww.stanford.edu/data/media/public/sep//prof/pvi.pdf Chapter 4.
	d = copy(in)
	n1,n2 = size(in)
	format = get(param,"format","angle") # output format: "angle" (in degrees wrt vertical) or "dip" (in y samples over x samples)
	pp = zeros(n1,n2); dx = wavekill(1.,pp,n1,n2,d)
	pp =  ones(n1,n2); dt = wavekill(0.,pp,n1,n2,d)
	dtdx = dt.*dx
	dxdx = dx.*dx
	dtdt = dt.*dt
	for i2 = 1 : n2
		dtdt[:,i2] = triangle(w1,n1,dtdt[:,i2])
		dxdx[:,i2] = triangle(w1,n1,dxdx[:,i2])
		dtdx[:,i2] = triangle(w1,n1,dtdx[:,i2])
	end
	for i1 = 1 : n1
		dtdt[i1,:] = triangle(w2,n2,dtdt[i1,:])
		dxdx[i1,:] = triangle(w2,n2,dxdx[i1,:])
		dtdx[i1,:] = triangle(w2,n2,dtdx[i1,:])
	end

	coh = sqrt((dtdx.*dtdx)./(dtdt.*dxdx))
	pp = -dtdx./dtdt

	coh = zeros(n1,n2); 
	for i1 = 1 : n1
		for i2 = 1 : n2
			if (abs(dtdt[i1,i2]) > 1e-6)
				coh[i1,i2] = sqrt((dtdx[i1,i2].*dtdx[i1,i2])./(dtdt[i1,i2].*dxdx[i1,i2]))
				pp[i1,i2] = -dtdx[i1,i2]./dtdt[i1,i2]
			else
				coh[i1,i2] = 0.
				pp[i1,i2]  = 0.			
			end
		end
	end
	for i2 = 1 : n2
		pp[:,i2] = triangle(w1,n1,pp[:,i2])
	end
	for i1 = 1 : n1
		pp[i1,:] = triangle(w2,n2,pp[i1,:])
	end

	res = wavekill(1.,pp,n1,n2,d)

	if (format == "angle")
		pp = atan(pp*dz/dx)*180/pi;
	end


	return coh,pp,res
end

function wavekill(aa,bb,n1,n2,uu)
	vv = zeros(n1,n2)
	s11 = -aa-bb; s12 = aa-bb; 
	s21 = -aa+bb; s22 = aa+bb;
	for i2 = 1 : n2-1
		for i1 = 1 : n1-1
			vv[i1,i2] = uu[i1  ,i2  ]*s11[i1,i2] + uu[i1  ,i2+1]*s12[i1,i2] + uu[i1+1,i2  ]*s21[i1,i2] + uu[i1+1,i2+1]*s22[i1,i2];
		end
	end
	vv[n1,:] = vv[n1-1,:]
	vv[:,n2] = vv[:,n2-1]	
	return vv
end

function triangle(nbox,nd,xx)
	yy = zeros(nd)
	pp = boxconv(nbox,nd,xx)
	np = nbox+nd-1
	qq = boxconv(nbox,np,pp)
	nq = nbox+np-1
	for i = 1 : nd
		yy[i] = qq[i+nbox-1]
	end
	for i = 1 : nbox - 1
		yy[i] = yy[i] + qq[nbox-i]
	end
	for i = 1 : nbox - 1
		yy[nd-i+1] = yy[nd-i+1] + qq[nd+(nbox-1)+i]
	end
	return yy
end

function boxconv(nbox,nx,xx)
	ny = nx+nbox-1
	yy = zeros(ny)
	bb = zeros(nx+nbox)
	bb[1] = xx[1]
	for i = 2 : nx
		bb[i] = bb[i-1] + xx[i]
	end
	for i = nx+1 : ny
		bb[i] = bb[i-1]
	end
	for i = 1 : nbox
		yy[i] = bb[i]
	end
	for i = nbox+1 : ny
		yy[i] = bb[i] - bb[i - nbox]
	end
	for i = 1 : ny
		yy[i] = yy[i]/nbox
	end
	return yy
end
