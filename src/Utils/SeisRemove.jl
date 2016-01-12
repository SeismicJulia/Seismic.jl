"""
**SeisRemove**

*Delete a seis file (deletes the text file, binary data file, and binary header file if there is one)*

**IN**
* in: input filename  

**OUT**  

*Credits: AS, 2015*

"""
function SeisRemove(filename::ASCIIString)

	filename_data = chomp(readall(`grep "in=" $filename` |> `tail -1` |> `awk '{print substr($1,5,length($1)-5) }' `))
	filename_headers = success(`grep "headers=" $filename`) ? chomp(readall(`grep "headers=" $filename` |> `tail -1` |> `awk '{print substr($1,10,length($1)-10) }' `)) : "NULL"
	rm(filename);
	rm(filename_data);
	if filename_headers != "NULL"
		rm(filename_headers);
	end	

end

function SeisRemove(filename::Array{ASCIIString,1})

	for j = 1 : length(in)
		SeisRemove(filename[j])
	end

end
