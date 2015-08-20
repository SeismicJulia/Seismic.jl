function ShotProfileLSWEM(m,d,param=Dict())

	# least squares shot profile wave equation migration of isotropic 1C data. 

	misfit = get(param,"misfit","misfit")   # misfit output text file
	wd = join(["tmp_LSM_wd_",string(int(rand()*100000))])
	CalculateSampling(d,wd)
	param["wd"] = wd
	param["operators"] = [ShotProfileWEM]
        ConjugateGradients(m,d,misfit,param)

end
