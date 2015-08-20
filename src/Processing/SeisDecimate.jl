function SeisDecimate(in,h,param)

	mode = get(param,"mode","random")
	style = get(param,"style","sxsygxgy")
	nt = size(in,1)
	nx = size(in,2)
	out = copy(in)
	
	if (mode=="random") # decimate data randomly
		perc = get(param,"perc",50)
		mask = rand(1,size(in,2)*size(in,3)*size(in,4)*size(in,5));
		mask[find(mask .< perc/100)] = 0;
		mask[find(mask .>= perc/100)] = 1;
		for it = 1 : nt
			out[it,:] = out[it,:].*mask;
		end
	else # decimate data regularly with respect to 4 spatial dimensions
		if (style == "sxsygxgy")
			key = ["t","isx","isy","igx","igy"]
		elseif (style=="mxmyhxhy")
			key = ["t","imx","imy","ihx","ihy"]
		elseif (style=="mxmyhaz")
			key = ["t","imx","imy","ih","iaz"]
		elseif (style=="sxsyhxhy")
			key = ["t","isx","isy","ihx","ihy"]
		elseif (style=="gxgyhxhy")
			key = ["t","igx","igy","ihx","ihy"]
		elseif (style=="sxsyhaz")
			key = ["t","isx","isy","ih","iaz"]
		elseif (style=="gxgyhaz")
			key = ["t","igx","igy","ih","iaz"]
		else
			error("style not defined.")
		end
		inc1 = get(param,"inc1",1)
		inc2 = get(param,"inc2",1)
		inc3 = get(param,"inc3",1)
		inc4 = get(param,"inc4",1)
		
		min_ix1 = getfield(h[1],symbol(key[2]))
		min_ix2 = getfield(h[1],symbol(key[3]))
		min_ix3 = getfield(h[1],symbol(key[4]))
		min_ix4 = getfield(h[1],symbol(key[5]))
		for ix = 1 : nx
			if ( getfield(h[ix],symbol(key[2]))/inc1 - abs(floor(getfield(h[ix],symbol(key[2]))/inc1)) != 0 )
				out[:,ix] = 0		
			end
			if ( getfield(h[ix],symbol(key[3]))/inc2 - abs(floor(getfield(h[ix],symbol(key[3]))/inc2)) != 0 )
				out[:,ix] = 0		
			end
			if ( getfield(h[ix],symbol(key[4]))/inc3 - abs(floor(getfield(h[ix],symbol(key[4]))/inc3)) != 0 )
				out[:,ix] = 0		
			end
			if ( getfield(h[ix],symbol(key[5]))/inc4 - abs(floor(getfield(h[ix],symbol(key[5]))/inc4)) != 0 )
				out[:,ix] = 0		
			end
		end
		
	end
	
    return out,h
end