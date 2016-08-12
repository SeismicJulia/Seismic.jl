function Pad5D(a,N1,N2,N3,N4,N5)

	n1,n2,n3,n4,n5 = size(a)
	b = zeros(Float32,N1,N2,N3,N4,N5)
	b[1:n1,1:n2,1:n3,1:n4,1:n5] = a

	return b
end
