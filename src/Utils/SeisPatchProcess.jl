@everywhere function MyProcess(filename)
	d1,h1,status = SeisRead(filename)
	for ifunc = 1 : length(f)
		func = f[ifunc]
		d2,h2 = func(d1,h1,f_param)
		d1 = copy(d2)
		h1 = copy(h2)
    end
	SeisWrite(filename,d1,h1)
	return(1)
end

@everywhere function SeisPatchProcess(in,out,param)	
# read data from disk, split into multidimensional overlapping patches, apply 
# processes, merge patches with tapers, and write to disk. Processing of patches 
# can be done in parallel by running the code using (for example): julia -p 2 script_name.jl
#
#
# Important Notice: You must make declare global variables f and f_param on every processor. You
# can do this in your main function by typing (for example): 
# @everywhere global f_param = ["style"=>"mxmyhxhy","Niter"=> 100,"alpha"=>1,"fmax"=>80]
# @everywhere global f = [SeisPOCS] 
#
# f is an array of functions that have the following syntax: d2,h2 = f(d1,h1,f_param), where 
# param is a dictionary (Dict) of parameters for the function.
#
# note that param should contain parameters for the patching and 
# unpatching operations.
#
# to execute the code type:  julia -p 4 my_code.jl where 4 can be replaced with the number of 
# processors you wish to use.
#
#
	patch_names = SeisPatch(in,out,param)
	a = pmap(MyProcess,patch_names)
	SeisUnPatch(patch_names,out,param)
	for ipatch = 1 : length(patch_names)
		rm(join([patch_names[ipatch] ".seisd" ]))
		rm(join([patch_names[ipatch] ".seish" ]))
	end

end
