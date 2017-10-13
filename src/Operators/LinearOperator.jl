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

function LinearOperator(m::AbstractString,d::AbstractString,operators,parameters;adj=true)

    rand_string = string(round(Int,rand()*100000))
    tmp_m = join(["tmp_CG_m_",rand_string])
    tmp_d = join(["tmp_CG_d_",rand_string])
    if adj
        SeisCopy(d,tmp_d)
        for j = 1 : 1 : length(operators)
            op = operators[j]
            op(tmp_m,tmp_d,true;parameters[j]...)
            SeisCopy(tmp_m,tmp_d)
        end
        SeisCopy(tmp_m,m)
        SeisRemove(tmp_m)
        SeisRemove(tmp_d)
    else
        SeisCopy(m,tmp_m)
        for j = length(operators) : -1 : 1
            op = operators[j]
            op(tmp_m,tmp_d,false;parameters[j]...)
            SeisCopy(tmp_d,tmp_m)
        end
        SeisCopy(tmp_d,d)
        SeisRemove(tmp_m)
        SeisRemove(tmp_d)
    end

end

function LinearOperator(m::Array{AbstractString,1},d::Array{AbstractString,1},operators,parameters;adj=true)

    rand_string = string(round(Int,rand()*100000))
    tmp_m = [join(["tmp_CG_m1_",rand_string]);join(["tmp_CG_m2_",rand_string])]
    tmp_d = [join(["tmp_CG_d1_",rand_string]);join(["tmp_CG_d2_",rand_string])]
    if adj
        SeisCopy(d,tmp_d)
        for j = 1 : 1 : length(operators)
            op = operators[j]
            if (length(methods(op,(Array,Array,Dict{Any,Any}))) > 0)
                op(tmp_m,tmp_d,true;parameters[j]...)
                SeisCopy(tmp_m,tmp_d)
            else
                for k = 1 : length(m)
                    op(tmp_m[k],tmp_d[k],true;parameters[j]...)
                end
                SeisCopy(tmp_m,tmp_d)
            end
        end
        SeisCopy(tmp_m,m)
        SeisRemove(tmp_m)
        SeisRemove(tmp_d)
    else
        SeisCopy(m,tmp_m)
        for j = length(operators) : -1 : 1
            op = operators[j]
            if (length(methods(op,(Array,Array,Dict{Any,Any}))) > 0)
                op(tmp_m,tmp_d,false;parameters[j]...)
                SeisCopy(tmp_d,tmp_m)
            else
                for k = 1 : length(m)
                    op(tmp_m[k],tmp_d[k],false;parameters[j]...)
                end
                SeisCopy(tmp_d,tmp_m)
            end
        end
        SeisCopy(tmp_d,d)
        SeisRemove(tmp_m)
        SeisRemove(tmp_d)
    end

end
