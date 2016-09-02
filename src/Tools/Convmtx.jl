"""
    Mc = convMtx(w, nt)

Generate convolution matrix

# Arguments
* `w :: Array{Float64,1}`: wavelet
* `nt :: Int64`          : trace length of seismic data

# Output
* `Mc :: Array{Float64,2}`: convolution matrix formulated by wavelet
"""
function convMtx(w::Array{Float64,1}, nr::Int64)
	  w  = vec(w)
	  nw = length(w)
		inds = floor(Int64, nw/2)
	  nd = nw + nr - 1
	  cm = zeros(nd, nr)
	  for i=1:nr
		    cm[i:i+nw-1, i] = w
	  end
	  cm = cm[inds+1:end-inds, :]
	  return cm
end

"""
    d = conv1(r, adj, w=Ricker())

generate synthetic data based on convolution model

# Arguments
* `r :: Array{Float64,2}`: reflectivities
* `adj :: Bool`          : switch between adjoint(true) and forward(false)

# Output
* `d :: Array{Float64,2}`: synthetic data
"""
function conv1(r::Array{Float64,2}, adj::Bool; w=Ricker())
    cm = convMtx(w, size(r,1))
		d  = zeros(r)
		if adj
			 cm = cm'
		end
		nt = size(r, 2)
		for it = 1 : nt
				d[:, it] = cm * r[:, it]
		end
	  return d
end
