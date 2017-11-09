function SeisSmoothGathers(m::AbstractString,d::AbstractString,adj;Nsmooth=3,Nrepeat=1)

	if (adj==false)
		m1,h,ext = SeisRead(m)
		ngather = size(m1,3)
		for igather = 1 : ngather
			A = m1[:,:,igather];
			A = smooth_angles(A;Nsmooth=Nsmooth,Nrepeat=Nrepeat)
			m1[:,:,igather] = A
		end
		SeisWrite(d,m1,h,ext)
	else
		d1,h,ext = SeisRead(d)
		ngather = size(d1,3)
		for igather = 1 : ngather
			A = d1[:,:,igather];
			A = smooth_angles(A;Nsmooth=Nsmooth,Nrepeat=Nrepeat)
			d1[:,:,igather] = A
		end
		SeisWrite(m,d1,h,ext)
	end

end

function smooth_angles(A;Nsmooth=3,Nrepeat=1)

	for iz = 1 : size(A,1)
		a = A[iz,:]
		for irepeat = 1 : Nrepeat
			a = mean_filter(a,Nsmooth)
		end
		A[iz,:] = a
	end

	return A;
end

function mean_filter(a,nw)

	b = vec(zeros(length(a),1))
	for ix = 1 : length(a)
		sum  = 0.
		nsum = 0.
		for iw = 1 : nw
			index1 = ix - Int(floor(nw/2)) - 1 + iw
			if (index1>0 && index1<length(a))
				sum += a[index1]
				nsum += 1.
			end
		end
		b[ix] = sum/nsum
	end
	return b
end
