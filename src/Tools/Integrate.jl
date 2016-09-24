function Integrate(nx, xx, yy, adj)

    t = 0.0f0
    if (adj)
	for i = nx : -1 : 1
	    t += yy[i]
	    xx[i] += t;
	end
    else
	for i = 1:nx
	    t += xx[i]
	    yy[i] += t;
	end
    end
    
end
