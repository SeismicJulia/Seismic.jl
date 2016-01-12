function LinearOperator(in,operators,parameters;adj=true)
	if adj
		d = copy(in)
		m = [];
		for j = 1 : 1 : length(operators)
			op = operators[j]
			m = op(d,true;parameters[j]...)
			d = copy(m)
		end
		return m
	else
		m = copy(in)
		d = [];
		for j = length(operators) : -1 : 1
			op = operators[j]
			d = op(m,false;parameters[j]...)
			m = copy(d)
		end
		return d
	end	

end


