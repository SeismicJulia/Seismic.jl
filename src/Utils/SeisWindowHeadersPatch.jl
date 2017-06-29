function SeisWindowHeadersPatch(in, out; key=[], minval=[], maxval=[], tmin=0,
                           tmax=99999, ntrace=1000000)



     itrace_in = 1
     itrace_out = 1
     NX = GetNumTraces(in)

     group="some"

     while itrace_in <= NX
         nx = NX - itrace_in + 1
	       ntrace = nx > ntrace ? ntrace : nx

	       h1 = SeisReadHeaders(in, group=group, key=key, itrace=itrace_in,
                                 ntrace=ntrace)
	       num_traces_in = length(h1)

         h2=WindowHeadersPatch(h1,key=key,minval=minval,maxval=maxval)
	       h1 = copy(h2)

         num_traces_out = length(h1)
	       SeisWriteHeaders(out, h1, itrace=itrace_out,
                             update_tracenum=false)
	       itrace_in += num_traces_in
	       itrace_out += num_traces_out
     end


    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
    filename_d_out = join([DATAPATH out "@data@"])
    filename_h_out = join([DATAPATH out "@headers@"])
    @compat nhead = length(fieldnames(Header))


    stream_h = open(filename_h_out)
    nx = round(Int,filesize(stream_h)/(nhead*4))

    h = GrabHeader(stream_h,1)
    close(stream_h)

    nt = convert(Int64,round((tmax - tmin)/h.d1)) + 1
    if nt > h.n1
	     nt = h.n1
    end

    minval = convert(Array{Float32,1},vec(minval))
    maxval = convert(Array{Float32,1},vec(maxval))
    key2 = String[]
    minval2 = Float32[]
    maxval2 = Float32[]
    for ikey=1:length(key)
      if key[ikey] != "t"
	       key2 = push!(key2,key[ikey])
	       minval2 = push!(minval2,minval[ikey])
	       maxval2 = push!(maxval2,maxval[ikey])
	    end
    end

    nx1=maxval2[1]-minval2[1]+1
    nx2=maxval2[2]-minval2[2]+1
    nx3=maxval2[3]-minval2[3]+1
    nx4=maxval2[4]-minval2[4]+1

    extent = Extent(convert(Int32,nt),
            convert(Int32,nx1), convert(Int32,nx2),
            convert(Int32,nx3), convert(Int32,nx4),
            convert(Float32,tmin),
            convert(Float32,1), convert(Float32,0), convert(Float32,0),
            convert(Float32,0), convert(Float32,h.d1),
            convert(Float32,1), convert(Float32,1), convert(Float32,1),
            convert(Float32,1), "Time", key[2], key[3], key[4], key[5], "s",
                    "index", "index", "index", "index", "")
    WriteTextHeader(out,extent,"native_float",4,filename_d_out,filename_h_out)

end

function WindowHeadersPatch(h_in;key=[],minval=[],maxval=[])

    minval = convert(Array{Float32,1},vec(minval))
    maxval = convert(Array{Float32,1},vec(maxval))
    nx = length(h_in)
    key2 = AbstractString[]
    minval2 = Float32[]
    maxval2 = Float32[]
    for ikey=1:length(key)
	     if key[ikey] != "t"
	        key2 = push!(key2,key[ikey])
	        minval2 = push!(minval2,minval[ikey])
	        maxval2 = push!(maxval2,maxval[ikey])
	     end
    end

    return RejectHeadersPatch(h_in,key2,minval2,maxval2,length(key2),length(h_in))

end

function RejectHeadersPatch(h_in::Array{Header,1}, key::Array{AbstractString,1},
                       minval::Array{Float32,1}, maxval::Array{Float32,1},
                       nkeys, nx)
    h_out = Header[]
    keep = true
    key_val = 0f0
    for j=1:nx
	     keep = true
	     for ikey=1:nkeys

	       key_val = convert(Float32,getfield(h_in[j],Symbol(key[ikey])))

	        if (key_val < minval[ikey] || key_val > maxval[ikey])
		          keep = false
	        end
	     end
	     if (keep==true)
	        h_out = push!(h_out,h_in[j])
	     end
    end
  return h_out
end
