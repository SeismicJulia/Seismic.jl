"""
**SeisCopy**

*Copy a seis file*

**IN**
* `in::AbstractString`: input filename
* `out::AbstractString`: output filename

**OUT**

*Credits: AS, 2015*

"""
function SeisCopy(in::AbstractString,out::AbstractString)

	filename_data_in = ParseDataName(in)
	filename_headers_in = ParseHeaderName(in)
	DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
	filename_data_out = join([DATAPATH out "@data@"])
	if filename_headers_in != "NULL"
		filename_headers_out = join([DATAPATH out "@headers@"])
	else
		filename_headers_out = "NULL"
	end
	extent = ReadTextHeader(in)
	WriteTextHeader(out,extent,"native_float",4,filename_data_out,filename_headers_out)
	run(`cp $filename_data_in $filename_data_out`)
	if filename_headers_in != "NULL"
		run(`cp $filename_headers_in $filename_headers_out`)
	end

end

function SeisCopy(in::Array{AbstractString,1},out::Array{AbstractString,1})

	for j = 1 : length(in)
		SeisCopy(in[j],out[j])
	end

end
