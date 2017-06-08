"""
    SeisSort(in, out;<keyword arguments>)

Sort a seis file using its header words

# Arguments
* `in`: input filename >> a text file with information about data extent, data and header file names; a binary file containing data and a binary file containing headers.
* `out`: output filename

# Keyword arguments
* `key=["imx","imy"]`
* `rev=false` : sort headers in decreasing order
* `ntrace=1000` : number of traces to read at a time

# Output
file `out` is created with data sorted.

*Credits: AS, 2015*
"""
function SeisSort(in, out;key=["imx","imy"],rev=false,ntrace=100000)
    filename_h = ParseHeaderName(in)
    stream_h = open(filename_h)
    seek(stream_h, header_count["n1"])
    nt = read(stream_h,Int32)
    nx = convert(Int64,filesize(stream_h)/(4*length(fieldnames(Header))))
    h = Header[]
    # find min and max for each key
    h1 = GrabHeader(stream_h,1)
    minval = Array(Float32,length(key))

    for ikey=1:length(key)
	minval[ikey] = getfield(h1,Symbol(key[ikey]))
    end
    for j=2:nx
	h1 = GrabHeader(stream_h,j)
	for ikey=1:length(key)
	    key_val = abs(getfield(h1,Symbol(key[ikey])))
	    if (key_val < minval[ikey])
		minval[ikey] = key_val
	    end
	end
    end
    mykey = vec(zeros(Float64,nx))
    seekstart(stream_h)
    for j=1:nx
	h1 = GrabHeader(stream_h,j)
	for ikey=1:length(key)
	    mykey[j] += ((getfield(h1,Symbol(key[ikey])) + minval[ikey] + 1)*
                         (10^(6*(length(key)-ikey))))
	end
    end
    close(stream_h)
    p = convert(Array{Int32,1},sortperm(mykey,rev=rev))
    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
    filename_d_out = join([DATAPATH out "@data@"])
    filename_h_out = join([DATAPATH out "@headers@"])
    nhead = length(fieldnames(Header))
    stream_h = open(filename_h_out)
    nx = Int(floor(filesize(stream_h)/(nhead*4)))
    close(stream_h)
    extent = ReadTextHeader(in)
    extent.n2 = nx
    extent.n3 = 1
    extent.n4 = 1
    extent.n5 = 1
    WriteTextHeader(out,extent,"native_float",4,filename_d_out,filename_h_out)
    FetchHeaders(filename_h,out,p,nx)
    Seismic.FetchTraces(in,out)
    tmp = join(["tmp_SeisSort_",string(Int(floor(rand()*100000)))])
    cp(out,tmp,remove_destination=true);
    @compat SeisProcessHeaders(out, tmp, [UpdateHeader],
                               [Dict(:itmin=>1,:itmax=>nt)])
    filename_h_tmp = join([DATAPATH tmp "@headers@"])
    filename_h_out = join([DATAPATH out "@headers@"])
    cp(filename_h_tmp,filename_h_out,remove_destination=true);
    rm(filename_h_tmp);
    rm(tmp);
end

function FetchHeaders(filename_h_in::AbstractString, filename_out::AbstractString,
                      p::Array{Int32,1}, nx)
    stream_h = open(filename_h_in)
    h = Header[]
    for j = 1:nx
	append!(h,[GrabHeader(stream_h,p[j])])
    end
    SeisWriteHeaders(filename_out,h,update_tracenum=false)
end
