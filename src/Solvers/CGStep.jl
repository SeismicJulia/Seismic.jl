function CGStep(m1,m2;a=1.,b=1.)
	# m1 = a*m1 + b*m2
	return a*m1 + b*m2
end

function CGStep(m1,m2,h1::Array{Header,1},h2::Array{Header,1};a=1.,b=1.)
	# m1 = a*m1 + b*m2
	return a*m1 + b*m2,h1
end

function CGStep(m1::AbstractString,m2::AbstractString;a=1.,b=1.)
	# m1 = a*m1 + b*m2

	#tmp = join(["tmp_CGStep_",string(round(Int,rand()*100000))])
	#@compat params = Dict(:a=>a,:b=>b)
	#SeisProcess(m1,m2,tmp,[CGStep],[params],key=["imx";"imy"])
	#SeisCopy(tmp,m1)
	#SeisRemove(tmp)
	d1,h,ext = SeisRead(m1)
	d2,h,ext = SeisRead(m2)
	SeisWrite(m1,a*d1[:,:] + b*d2[:,:],h,ext)
end

function CGStep(m1::Array{AbstractString,1},m2::Array{AbstractString,1};a=[1. 1. 1.],b=[1. 1. 1.])
	# m1 = a*m1 + b*m2
	for j = 1 : length(m1)
		CGStep(m1[j],m2[j],a=a[j],b=b[j])
	end

end
