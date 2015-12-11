function rectangle_filter(nb,nx,xx,yy)

	bb = zeros(Float32,nx+nb)
	ny = nx+nb-1
	bb[1] = xx[1]
	for i = 2 : nx
		bb[i] = bb[i-1] + xx[i]
	end
	for i = nx+1 : ny
		bb[i] = bb[i-1]
	end
	for i = 1 : nb
		yy[i] = bb[i]
	end
	for i = nb+1 : ny
		yy[i] = bb[i] - bb[i-nb]
	end
	for i = 1 : ny
		yy[i] = yy[i] / nb
	end

end

function triangle_filter(nr, m1, n12, uu, vv)

	pp = zeros(Float32,n12+nr-1)
	qq = zeros(Float32,n12+nr+nr-2)
	tt = zeros(Float32,n12)
	for i = 1 : n12
		qq[i] = uu[i]
	end
	rectangle_filter(nr,n12,qq,pp)
	rectangle_filter(nr,nr+n12-1,pp,qq)
	for i = 1 : n12
		tt[i] = qq[i+nr-1]
	end
	for i = 1 : nr-1
		tt[i] = tt[i] + qq[nr-i]
	end
	for i = 1 : nr-1
		tt[n12-i+1] = tt[n12-i+1] + qq[n12+(nr-1)+i]
	end
	for i = 1 : n12
		vv[i] = tt[i]
	end

end

function integrate(nx,xx,yy,adj)

	t = 0.0f0
	if (adj)
		for i = nx : -1 : 1
			t += yy[i]
			xx[i] += t;
		end
	else
		for i = 1 : nx
			t += xx[i]
			yy[i] += t;
		end
	end
end
