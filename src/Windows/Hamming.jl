"""
    Hamming(n)

Create a Hamming window.

# Arguments
* `n::Int`: number of samples

# Example
```julia
julia> w = Hamming(51); plot(w);
```
"""

function Hamming(n::Int)

    n==1 && return [1.]
    n -= 1
    k = collect(0:1:n)
    hamm = 0.54 - 0.46*cos.(2pi.*k/n)

end
