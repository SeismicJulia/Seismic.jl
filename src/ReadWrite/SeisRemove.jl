"""
**SeisRemove**

*Delete a seis file (deletes the text file, binary data file, and binary header file if there is one)*

**IN**
* in: input filename  

**OUT**  

*Credits: AS, 2015*

"""
function SeisRemove(filename::ASCIIString)

	filename_data = ParseDataName(filename)
	filename_headers = ParseHeaderName(filename)	
	rm(filename);
	rm(filename_data);
	if filename_headers != "NULL"
		rm(filename_headers);
	end	

end

function SeisRemove(filename::Array{ASCIIString,1})

	for j = 1 : length(filename)
		SeisRemove(filename[j])
	end

end
