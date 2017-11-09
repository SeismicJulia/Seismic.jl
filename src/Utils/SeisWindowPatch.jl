"""
**SeisWindowPatch**
*Window a seis file using header words.*
**IN**
* in
* out
* key=[]
* minval=[]
* maxval=[]
* note that windowing along the time axis is achieved by using the key "t".
**OUT**
*Credits: AS, FC, 2017*
"""

function SeisWindowPatch(in::AbstractString,out::AbstractString;key=[],minval=[],maxval=[],it_nt=9e9)
  println("processing patch ",out)
  DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
  extent = ReadTextHeader(in)
  tmin = extent.o1
  tmax = extent.o1 + extent.d1*(extent.n1-1)
  for ikey=1:length(key)
	   if key[ikey] == "t"
	      tmin = minval[ikey]
	      tmax = maxval[ikey]
	   end
  end

  #itmin is the sample at which tmin is saved
  itmin = convert(Int32,round((tmin - extent.o1)/extent.d1) + 1)

  if (itmin < 1)
	   itmin = 1
  end
  itmax = convert(Int32,round((tmax - extent.o1)/extent.d1) + 1)
  if (itmax) > extent.n1
	   itmax = extent.n1
  end


  SeisWindowHeadersPatch(in, out; key=key, minval=minval, maxval=maxval, tmin=tmin,
                      tmax=tmax,it_nt=it_nt)
  FetchTracesPatch(in,out;itmin=itmin,itmax=itmax ,key=key, minval=minval, maxval=maxval)
  tmp = join(["tmp_SeisWindow_",string(round(Int,rand()*100000))])

 #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 #@compat SeisProcessHeaders(out, tmp, [UpdateHeader],[Dict(:itmin=>itmin,:itmax=>itmax)])

  itrace_in = 1
  itrace_out = 1
  ntrace=1000000
  NX = GetNumTraces(out)
  while itrace_in <= NX
    nx = NX - itrace_in + 1
	  ntrace = nx > ntrace ? ntrace : nx
	  h1 = SeisReadHeaders(out, group="some", itrace=itrace_in,
                                 ntrace=ntrace)
    ext=ReadTextHeader(out)

	  num_traces_in = length(h1)
	  h2=UpdateHeaderPatch(h1;itrace=itrace_in,itmin=itmin,itmax=itmax)
	  h1 = copy(h2)
	  num_traces_out = length(h1)
	  SeisWriteHeaders(tmp, h1, itrace=itrace_out)
    itrace_in += num_traces_in
	  itrace_out += num_traces_out
  end

  filename_h_tmp = join([DATAPATH tmp "@headers@"])
  filename_h_out = join([DATAPATH out "@headers@"])
  cp(filename_h_tmp,filename_h_out,remove_destination=true);
  rm(filename_h_tmp);
    #rm(tmp);
end


########################################################################################
#####################################################################################
function FetchTracesPatch(in::AbstractString, out::AbstractString; ntrace=500, itmin=round(Int,1),
    itmax=round(Int,9999999999),key=[], minval=[], maxval=[])

    NX = GetNumTraces(out)

    itrace = round(Int,1)
    filename_data_in = ParseDataName(in)

    DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
    filename_data_out = join([DATAPATH out "@data@"])

    stream_in  = open(filename_data_in)

    stream_out  = open(filename_data_out,"w")
    extent = ReadTextHeader(in)
    nt = extent.n1

    if itmin < 1
      itmin = round(Int,1)
    end
    if itmax > nt
      itmax = round(Int,nt)
    end

    @compat nhead = length(fieldnames(Header))

    while itrace <= NX
      if (itrace > 1)
	       stream_out = open(filename_data_out,"a")
	    end
      nx = NX - itrace + 1

      ntrace = nx > ntrace ? ntrace : nx
      h  = SeisReadHeaders(out,group="some",itrace=itrace,ntrace=ntrace)
      d = zeros(Float32,itmax-itmin+1,ntrace)

      SeekTracesPatch!(d,stream_in,h,itmin,itmax,nt,round(Int,ntrace))


    write(stream_out,d)


    close(stream_out)
    itrace += ntrace
  end
  close(stream_in)
end

###############################################################################

function SeekTracesPatch!{T}(d::AbstractArray{T,2}, stream_in::IOStream,
                        h::Array{Header,1},itmin,itmax,nt,ntrace)

    d1 = zeros(Float32,nt)
    for ix = 1 : ntrace
	     position = 4*nt*(h[ix].tracenum-1)
	     seek(stream_in,position)
	     d1 = read(stream_in,Float32,nt)
	     for it = itmin : itmax
	        d[it - itmin + 1,ix] = d1[it]
	     end
    end
    nothing
end


##############################################################################
#function SeekTracesPatch{T}(d::AbstractArray{T,2}, h::Array{Header,1},out::String, stream_in::IOStream,stream_out::IOStream, itmin,itmax,nt,ntrace;key=[],minval=[],maxval=[])
# println(itmin," ",itmax," ",nt)
#    d1 = zeros(Float32,nt)
#extent_out=ReadTextHeader(out)
#nx2=extent_out.n3
#nx3=extent_out.n4
#nx4=extent_out.n5
#    for ix = 1 : ntrace
#	position = 4*nt*(h[ix].tracenum-1)
#	seek(stream_in,Int(position))

#	d1 = read(stream_in,Float32,nt)
#	for it = itmin : itmax
#	    d[it - itmin + 1,ix] = d1[it]

#	end



#        itrace_out=(h[ix].imx-minval[2])*nx2*nx3*nx4+(h[ix].imy-minval[3])*nx3*nx4+(h[ix].ihx-minval[4])*nx4+h[ix].ihy-minval[5]+1

#        position_out=4*nt*(itrace_out-1)
#println(out," ",h[ix].imx," ",minval[2]," ",itrace_out)
#        seek(stream_out,Int(position_out))
#        write(stream_out,convert(Array{Float32,1},d[1:itmax-itmin+1,ix]))
#    end
#    nothing
#end

##############################################################################
##############################################################################
function UpdateHeaderPatch(h;itrace=1,itmin=1,itmax=9e9)
    ot = (itmin-1)*h[1].d1 + h[1].o1

    nt = itmax - itmin + 1
    for j = 1 : length(h)
	h[j].o1 = ot
	h[j].n1 = nt
	h[j].tracenum = j + itrace - 1
    end
    return h
end
