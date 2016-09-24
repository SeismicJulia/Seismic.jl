function SeisDecimate(in;mode="random",perc=50,incx1=1,incx2=1,incx3=1,incx4=1)
    
    if (mode=="random") # decimate data randomly
        out = copy(in)
	mask = rand(1,size(in,2)*size(in,3)*size(in,4)*size(in,5));
	mask[find(mask .< perc/100)] = 0;
	mask[find(mask .>= perc/100)] = 1;
	for it = 1 : size(in,1)
	    out[it,:] = out[it:it,:].*mask;
	end
    else # decimate data regularly with respect to 4 spatial dimensions
        out = zeros(in)
	if (size(in),5) > 1
	    out[:,1:incx1:end,1:incx2:end,1:incx3:end,1:incx4:end] = in[:,1:incx1:end,1:incx2:end,1:incx3:end,1:incx4:end]
	elseif (size(in),4) > 1
            out[:,1:incx1:end,1:incx2:end,1:incx3:end] = in[:,1:incx1:end,1:incx2:end,1:incx3:end]
        elseif (size(in),3) > 1
            out[:,1:incx1:end,1:incx2:end] = in[:,1:incx1:end,1:incx2:end]
        elseif (size(in),2) > 1
            out[:,1:incx1:end] = in[:,1:incx1:end]
        end		
	
    end
    
    return out
end
