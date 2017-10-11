function TriangleFilter(nr, n12, uu, vv)

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
