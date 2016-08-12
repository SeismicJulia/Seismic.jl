function Convmtx(w,nr)
	w  = vec(w)
	nw = length(w)
	nd = nw + nr - 1
	cm = zeros(nd, nr)
	for i=1:nr
		cm[i:i+nw-1, i] = w
	end
	cm = sparse(cm)
	return cm
end

