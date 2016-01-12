"""
**SeisCopy**

*Copy a seis file*

**IN**
* in: input filename  
* out: output filename  

**OUT**  

*Credits: AS, 2015*

"""
function SeisCopy(in::ASCIIString,out::ASCIIString)
	filename_data_in = chomp(readall(`grep "in=" $in` |> `tail -1` |> `awk '{print substr($1,5,length($1)-5) }' `))
	filename_headers_in = success(`grep "headers=" $in`) ? chomp(readall(`grep "headers=" $filename` |> `tail -1` |> `awk '{print substr($1,10,length($1)-10) }' `)) : "NULL"
	DATAPATH = get(ENV,"DATAPATH","./")
	filename_data_out = join([DATAPATH out "@data@"])
	if filename_headers != "NULL"	
		filename_headers_out = join([DATAPATH out "@headers@"])
	else
		filename_headers_out = "NULL"
	end		
	extent = ReadTextHeader(in)
	esize = int(chomp(readall(`grep "esize" $filename` |> `tail -1` |> `awk '{print substr($1,7,length($1))}'` )))
	WriteTextHeader(out,extent,"native_float",esize,filename_data_out,filename_headers_out)
	run(`cp $filename_data_in $filename_data_out`)
	if filename_headers != "NULL"
		run(`cp $filename_headers_in $filename_headers_out`)
	end	

end

function SeisCopy(in::Array{ASCIIString,1},out::Array{ASCIIString,1})

	for j = 1 : length(in)
		SeisCopy(in[j],out[j])
	end
	
end
