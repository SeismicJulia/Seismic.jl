"""
    dpp = SeisConv(vp, dx, dt, f0)

Generate synthetic data based on convolution model

# Arguments
* `vp :: Array{Float64,2}`: 2D p-wave velocity model
* `dx :: Float64`       : sptial interval of velocty model
* `dt :: Float64`       : time sample
* `f0 :: Float64`       : dominant frequency of Ricker wavelet

# Output
* `dpp :: Array{Float64,2}`: synthetic data based on convolution model
"""
function SeisConv(vp::Array{Float64,2}, dx::Float64, dt::Float64, f0::Float64)
    (m, n) = size(vp)
    nr = m-1
    rc_pp = zeros(nr, n)
    for ix = 1:n
        for k = 1:nr
            rc_pp[k, ix] = (vp[k+1,ix]-vp[k,ix]) / (vp[k+1,ix]+vp[k,ix])
        end
    end
    tpp = zeros(nr, n)
#   the time location of reflecter
    for ix = 1:n
        for k = 1:nr
            for l =1: k
                tpp[k, ix] = tpp[k, ix] + 2*dx/vp[l, ix]
            end
        end
    end
    ntp = floor(Int64, maximum(tpp)/dt) + 30
		rpp = zeros(ntp, n)
    for ix = 1 : n
        for k = 1 : nr
            jp = floor(Int64, tpp[k,ix]/dt) + 1
            if abs(rpp[jp, ix]) < abs(rc_pp[k, ix])
               rpp[jp,ix] = rc_pp[k, ix]
            end
        end
    end
    w = Ricker(f0, dt)
		dpp = conv1(rpp, false, w=w)
    return dpp
end
