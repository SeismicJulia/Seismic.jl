"""
    PadFirstAxis(d, N1)

Zero-padding first axis of N-dimensional array `d`.

# Arguments
* `d::Array{Real, N}`: N-dimensional data.
* `N1::Int`: total number of samples for the first axis of the zero-padded data.

# Examples
```julia
julia> w = Ricker(); wpad = PadFirstAxis(w, 128); plot(w); plot(wpad+1.0)

julia> d, ext = SeisHypEvents(); dpad = PadFirstAxis(d, 512);
SeisPlot(d); SeisPlot(dpad)
```
"""

function PadFirstAxis{T<:Real, N}(a::Array{T,N}, N1::Int)
    dims = size(a)
    n1 = dims[1]
    N1 > n1 || error("N1 must be greater than size(a,1)")
    length(dims) == 1 ? nx = 1 : nx = *(dims[2:end]...)
    b = zeros(T, N1, nx)
    for j = 1:nx
	b[1:n1, j] = a[:, j]
    end
    reshape(b, N1, dims[2:end]...)
end
