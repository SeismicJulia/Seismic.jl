include("Header.jl")

function SeisHeaderInfo(in)
#
#   print Seis header information
#

    key = ["tracenum","o1","n1","d1","sx","sy","gx","gy","mx","my","hx","hy","h","az","ang","isx","isy","igx","igy","imx","imy","ihx","ihy","ih","iaz","iang","selev","gelev","sstat","gstat","trid"]

    filename = join([in ".seish"])
    stream = open(filename)
    nhead = 27
    nx = int(filesize(stream)/124)
    h = GrabHeader(stream,1)
    println("Displaying information for ", filename," (",nx," traces):")
    min_h = vec(zeros(Float32,length(key)))
    max_h = vec(zeros(Float32,length(key)))
    mean_h = vec(zeros(Float32,length(key)))
    for ikey=1:length(key)
    	min_h[ikey] = getfield(h,symbol(key[ikey]))
    	max_h[ikey] = getfield(h,symbol(key[ikey]))
    	mean_h[ikey] += getfield(h,symbol(key[ikey]))
    end
    for j=2:nx
    	#println("j=",j," nx=",nx)
    	h = GrabHeader(stream,j)
        for ikey=1:length(key)
        	key_val = getfield(h,symbol(key[ikey]))
            if (key_val < min_h[ikey])
        		min_h[ikey] = key_val
            end
            if (key_val > max_h[ikey])
        		max_h[ikey] = key_val
            end
            mean_h[ikey] += key_val
        end
    end
    for ikey=1:length(key)
    	mean_h[ikey] /= nx
    end
    close(stream)
    println("       Key          Minimum          Maximum             Mean");    
    println("=============================================================")
    for ikey=1:length(key)
		@printf("%10s      %11.3f      %11.3f      %11.3f\n",key[ikey],min_h[ikey],max_h[ikey],mean_h[ikey])
    end
    
end
