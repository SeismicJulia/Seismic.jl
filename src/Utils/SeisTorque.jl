function bind_torque_procs()

	home = ENV["HOME"]

	# find out what cluster we are on
	node_file_name = ENV["PBS_NODEFILE"]

	# read PBS_NODEFILE
	filestream = open(node_file_name)
	seekstart(filestream)
	node_file = readlines(filestream)

	# strip eol
	node_file = map(x->strip(x,['\n']),node_file)

	# get number of workers on each node
	procs = Dict{ASCIIString,Int}()
	for n in node_file
		procs[n] = get(procs,n,0) + 1
	end

	println("name of compute nodes and number of workers:")
	println(procs)

	# add processes on master itself
	master = ENV["HOSTNAME"]

	if procs[master] > 1
		addprocs(procs[master]-1)
		println("added $(procs[master]-1) processes on master itself")
	end

	# get a machine file for other hosts
	machines = ASCIIString[]
	for i in 1:length(node_file)
		if node_file[i] != master
			push!(machines,node_file[i])
		end
	end


	println("processes on worker machines:")
	println(machines)

	# add to julia home
	println("adding machines to JULIA_HOME: $JULIA_HOME")
	for m in machines
		addprocs([m], dir= JULIA_HOME)
	end

	println("done")
end

function bind_torque_procs(ppn::Int)

	# add only ppn processes per node

	home = ENV["HOME"]

	# find out what cluster we are on
	node_file_name = ENV["PBS_NODEFILE"]

	# read PBS_NODEFILE
	filestream = open(node_file_name)
	seekstart(filestream)
	node_file = readlines(filestream)

	# strip eol
	node_file = map(x->strip(x,['\n']),node_file)

	# get number of workers on each node
	procs = Dict{ASCIIString,Int}()
	for n in node_file
		procs[n] = get(procs,n,0) + 1
	end

	println("name of compute nodes and number of workers:")
	println(procs)

	# add processes on master itself
	master = ENV["HOSTNAME"]

	wrker = 0
	while wrker < ppn
		addprocs(1)
		wrker += 1
	end
	println("added $wrker processes on master itself")

	# add procs on other machines
	for (k,v) in procs
		wrker = 0
		if k!=master
			while wrker < ppn
				addprocs([k], dir= JULIA_HOME)
				# println("addprocs($k)")
				wrker += 1
			end
			println("added $wrker processes on machine $k")
		end
	end

	println("done")
end

