function PadFirstAxis{T<:Real, N}(a::Array{T,N}, N1::Int)
    dims = size(a)
    n1 = dims[1]
    length(dims) == 1 ? nx = 1 : nx = *(dims[2:end]...)
    b = zeros(T, N1, nx)
    for j = 1:nx
	b[1:n1, j] = a[:, j]
    end
    reshape(b, N1, dims[2:end]...)
end
