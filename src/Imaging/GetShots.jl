function GetShots(in)

	filename_h = join([filename ".seish"])
	stream = open(filename_h)    
	sx = Array(Any)
	sy = Array(Any)
	key = ["sx","sy"]
	curr = zeros(length(key),1)
	prev = copy(curr)
	ntraces = int(filesize(stream)/(4*length(names(Header))))
	nsx = 1
	h1 = GrabHeader(stream,1)
	for ikey=1:length(key)
		curr[ikey] = getfield(h1,Symbol(key[ikey]))
	end
	prev = copy(curr)
	sx = push!(sx,h1.sx)
	sy = push!(sy,h1.sy)
	for j=2:ntraces
		h1 = GrabHeader(stream,j)
		for ikey=1:length(key)
			curr[ikey] = getfield(h1,Symbol(key[ikey]))
		end
		if curr != prev
			nx = j - itrace
			sx = push!(sx,h1.sx)
			sy = push!(sy,h1.sy)
		end
		prev = copy(curr)
	end
	close(stream)
	return sx,sy

end
