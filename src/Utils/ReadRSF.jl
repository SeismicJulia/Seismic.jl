using Seismic

function ReadRSF(filename;swap_bytes=false)

	rsf_filename = chomp(readall(`grep "in" $filename` |> `grep "@"` |> `tail -1` |> `awk '{print substr($1,5,length($1)-5) }' `))
	n1 = success(`grep "n1" $filename`) ? int(chomp(readall(`grep "n1" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	n2 = success(`grep "n2" $filename`) ? int(chomp(readall(`grep "n2" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	n3 = success(`grep "n3" $filename`) ? int(chomp(readall(`grep "n3" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	n4 = success(`grep "n4" $filename`) ? int(chomp(readall(`grep "n4" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	n5 = success(`grep "n5" $filename`) ? int(chomp(readall(`grep "n5" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	o1 = success(`grep "o1" $filename`) ? float(chomp(readall(`grep "o1" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	o2 = success(`grep "o2" $filename`) ? float(chomp(readall(`grep "o2" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	o3 = success(`grep "o3" $filename`) ? float(chomp(readall(`grep "o3" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	o4 = success(`grep "o4" $filename`) ? float(chomp(readall(`grep "o4" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	o5 = success(`grep "o5" $filename`) ? float(chomp(readall(`grep "o5" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 0
	d1 = success(`grep "d1" $filename`) ? float(chomp(readall(`grep "d1" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	d2 = success(`grep "d2" $filename`) ? float(chomp(readall(`grep "d2" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	d3 = success(`grep "d3" $filename`) ? float(chomp(readall(`grep "d3" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	d4 = success(`grep "d4" $filename`) ? float(chomp(readall(`grep "d4" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	d5 = success(`grep "d5" $filename`) ? float(chomp(readall(`grep "d5" $filename` |> `tail -1` |> `awk '{print substr($1,4,length($1))}'` ))) : 1
	label1 = success(`grep "label1" $filename`) ? chomp(readall(`grep "label1" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	label2 = success(`grep "label2" $filename`) ? chomp(readall(`grep "label2" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	label3 = success(`grep "label3" $filename`) ? chomp(readall(`grep "label3" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	label4 = success(`grep "label4" $filename`) ? chomp(readall(`grep "label4" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	label5 = success(`grep "label5" $filename`) ? chomp(readall(`grep "label5" $filename` |> `tail -1` |> `awk '{print substr($0,10,length($0)-10)}'` )) : ""
	unit1 = success(`grep "unit1" $filename`) ? chomp(readall(`grep "unit1" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	unit2 = success(`grep "unit2" $filename`) ? chomp(readall(`grep "unit2" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	unit3 = success(`grep "unit3" $filename`) ? chomp(readall(`grep "unit3" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	unit4 = success(`grep "unit4" $filename`) ? chomp(readall(`grep "unit4" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	unit5 = success(`grep "unit5" $filename`) ? chomp(readall(`grep "unit5" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	title = success(`grep "title" $filename`) ? chomp(readall(`grep "title" $filename` |> `tail -1` |> `awk '{print substr($0,9,length($0)-9)}'` )) : ""
	stream = open(rsf_filename)
	dtype = chomp(readall(`grep "data_format" $filename` |> `tail -1` |> `awk '{print substr($1,14,length($1)-14)}'` ))
	dtype = dtype == "native_float" ? Float32 : Complex{Float32}
	esize = int(chomp(readall(`grep "esize" $filename` |> `tail -1` |> `awk '{print substr($1,7,length($1))}'` )))
	total = int(filesize(stream)/esize)
	d = read(stream,dtype,total)
	#if (swap_bytes==true)
	#	d = bswap_vector(d)
	#end
	close(stream)
	extent = Extent(convert(Int32,n1),convert(Int32,n2),convert(Int32,n3),convert(Int32,n4),convert(Int32,n5),
		   convert(Float32,o1),convert(Float32,o2),convert(Float32,o3),convert(Float32,o4),convert(Float32,o5),
		   convert(Float32,d1),convert(Float32,d2),convert(Float32,d3),convert(Float32,d4),convert(Float32,d5),
		   label1,label2,label3,label4,label5,
		   unit1,unit2,unit3,unit4,unit5,
		   title)
	if n5 == 1 && n4 == 1 && n3 == 1 && n2 == 1 
		d = reshape(d,n1)
	elseif n5 == 1 && n4 == 1 && n3 == 1
		d = reshape(d,n1,n2)
	elseif n5 == 1 && n4 == 1
		d = reshape(d,n1,n2,n3)
	elseif n5 == 1
		d = reshape(d,n1,n2,n3,n4)
	else
		d = reshape(d,n1,n2,n3,n4,n5)
	end
	
	return d,extent
end

