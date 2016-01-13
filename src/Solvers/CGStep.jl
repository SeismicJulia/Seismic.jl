function CGStep(m1,m2;a=1.,b=1.)
	# m1 = a*m1 + b*m2
	return a*m1 + b*m2,h1
end

function CGStep(m1::ASCIIString,m2::ASCIIString;a=1.,b=1.)
	# m1 = a*m1 + b*m2	

	tmp = join(["tmp_CGStep_",string(int(rand()*100000))])
	SeisProcess(m1,m2,tmp;group="some",f="[CGStep]",ntrace=100000)
	SeisCopy(tmp,m1)
	SeisRemove(tmp)
end

function CGStep(m1::Array{ASCIIString,1},m2::Array{ASCIIString,1};a=[1. 1. 1.],b=[1. 1. 1.])
	# m1 = a*m1 + b*m2	
	for j = 1 : length(m1)
		CGStep(m1[j],m2[j],a=a[j],b=b[j])
	end
	
end

