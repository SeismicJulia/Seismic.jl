
<a id='Utils-1'></a>

# Utils


Seismic.jl provides tools to work with real data sets like binning, patching, unpatching, data sorting and windowing. 


<a id='SeisGeometry-1'></a>

## SeisGeometry

<a id='Seismic.SeisGeometry' href='#Seismic.SeisGeometry'>#</a>
**`Seismic.SeisGeometry`** &mdash; *Function*.



```
SeisGeometry(in;<keyword arguments>)
```

Update headers with geometry information.

**Arguments**

  * `in`: input filename

**Keyword arguments**

  * `ang=90`: inline direction measured in degrees CC from East
  * `gamma=1`: vp/vs ratio for PS Asymptotic Conversion Point gathers (use gamma=1 for PP data)
  * `osx=0`,`osy=0`,`ogx=0`,`ogy=0` : origin for source and receiver coordinate system
  * `omx=0`,`omy=0`,`ohx=0`,`ohy=0`: origin for midpoint and offset coordinate system
  * `oaz=0`,`oh=0` : origin for azimuth and offset coordinate system
  * `dsx=1`,`dsy=1`,`dgx=1`,`dgy=1`: source and receiver step-size
  * `dmx=1`,`dmy=1`,`dhx=1`,`dhy=1`: midpoint and offset step-size
  * `dh=1`,`daz=1`: offset and azimuth step-size

**Outputs**

the .seish file is updated with the following information:

  * hx,hy,h,az,mx,my : calculated offset, azimuth and midpoint
  * isx,isy,igx,igy,imx,imy,ihx,ihy,ih,iaz: calculated grid nodes for source and receiver position and midpoint, offset and azimuth.

*Credits: A. Stanton, F.Carozzi,2017*


<a target='_blank' href='https://github.com/fercarozzi/myseismicjulia/tree/c2832f8331d8b4cba573c54c2dd0183c518801d7/src/Utils/SeisGeometry.jl#L1-L27' class='documenter-source'>source</a><br>


##SeisBinHeaders

<a id='Seismic.SeisBinHeaders' href='#Seismic.SeisBinHeaders'>#</a>
**`Seismic.SeisBinHeaders`** &mdash; *Function*.



```
SeisBinHeaders(in,out; <keyword arguments>)
```

Sequentially bin seismic headers using the available grid information.

Keyword arguments should be consistent with SeisGeometry keyword arguments.

**Arguments**

  * `in`: filename of input, irregularly sampled data
  * `out`: filename of output, regularly sampled data

**Keyword arguments**

  * `style="sxsygxgy"`: bin style. Options: "mxmyhxhy","mxmyhaz","sxsyhxhy","gxgyhxhy","sxsyhaz","gxgyhaz"
  * `ang=90`: inline direction measured in degrees CC from East
  * `gamma=1`: vp/vs ratio for PS Asymptotic Conversion Point gathers (use gamma=1 for PP data)
  * `osx=0`,`osy=0`,`ogx=0`,`ogy=0` : origin for source and receiver coordinate system
  * `omx=0`,`omy=0`,`ohx=0`,`ohy=0`: origin for midpoint and offset coordinate system
  * `oaz=0`,`oh=0` : origin for azimuth and offset coordinate system
  * `dsx=1`,`dsy=1`,`dgx=1`,`dgy=1`: source and receiver step-size
  * `dmx=1`,`dmy=1`,`dhx=1`,`dhy=1`: midpoint and offset step-size
  * `dh=1`,`daz=1`: offset and azimuth step-size
  * `min_isx=0`,`max_isx=0`,`min_isy=0`,`max_isy=0`: grid extreme values for sources
  * `min_igx=0`,`max_igx=0`,`min_igy=0`,`max_igy=0`: grid extreme values for receivers
  * `min_imx=0`,`max_imx=0`,`min_imy=0`,`max_imy=0`: grid extreme values for midpoints
  * `min_ihx=0`,`max_ihx=0`,`min_ihy=0`,`max_ihy=0`: grid extreme values for offsets
  * `min_ih=0`,`max_ih=0`,`min_iaz=0`,`max_iaz=0`: grid extreme values for azimuth and offset
  * `ntrace=10000`: maximum number of traces processed at a time

**Output**

In file `out`, binned headers are created.

*Credits: Aaron Stanton, F Carozzi,2017*


<a target='_blank' href='https://github.com/fercarozzi/myseismicjulia/tree/c2832f8331d8b4cba573c54c2dd0183c518801d7/src/Utils/SeisBinHeaders.jl#L1-L34' class='documenter-source'>source</a><br>


##SeisBinData

<a id='Seismic.SeisBinData' href='#Seismic.SeisBinData'>#</a>
**`Seismic.SeisBinData`** &mdash; *Function*.



```
SeisBinData(in,out; <keyword arguments>)
```

Sequentially bin seismic data using already binned trace headers (SeisBinHeaders). Input arguments should be consistent with SeisBinHeaders input arguments.

**Arguments**

  * `in`: filename of input, irregularly sampled data
  * `out`: filename of output, regularly sampled data

**Keyword arguments**

  * `style="sxsygxgy"`: bin style. Options: "mxmyhxhy","mxmyhaz","sxsyhxhy","gxgyhxhy","sxsyhaz","gxgyhaz"
  * `ang=90`: inline direction measured in degrees CC from East
  * `gamma=1`: vp/vs ratio for PS Asymptotic Conversion Point gathers (use gamma=1 for PP data)
  * `osx=0`,`osy=0`,`ogx=0`,`ogy=0` : origin for source and receiver coordinate system
  * `omx=0`,`omy=0`,`ohx=0`,`ohy=0`: origin for midpoint and offset coordinate system
  * `oaz=0`,`oh=0` : origin for azimuth and offset coordinate system
  * `dsx=1`,`dsy=1`,`dgx=1`,`dgy=1`: source and receiver step-size
  * `dmx=1`,`dmy=1`,`dhx=1`,`dhy=1`: midpoint and offset step-size
  * `dh=1`,`daz=1`: offset and azimuth step-size
  * `min_isx=0`,`max_isx=0`,`min_isy=0`,`max_isy=0`: grid extreme values for sources
  * `min_igx=0`,`max_igx=0`,`min_igy=0`,`max_igy=0`: grid extreme values for receivers
  * `min_imx=0`,`max_imx=0`,`min_imy=0`,`max_imy=0`: grid extreme values for midpoints
  * `min_ihx=0`,`max_ihx=0`,`min_ihy=0`,`max_ihy=0`: grid extreme values for offsets
  * `min_ih=0`,`max_ih=0`,`min_iaz=0`,`max_iaz=0`: grid extreme values for azimuth and offset
  * `ntrace=10000`: maximum number of traces processed at a time

**Output**

In file `out`, the binned data is created.

*Credits: Aaron Stanton, Fernanda Carozzi, 2017*


<a target='_blank' href='https://github.com/fercarozzi/myseismicjulia/tree/c2832f8331d8b4cba573c54c2dd0183c518801d7/src/Utils/SeisBinData.jl#L1-L33' class='documenter-source'>source</a><br>


##SeisPatch

<a id='Seismic.SeisPatch' href='#Seismic.SeisPatch'>#</a>
**`Seismic.SeisPatch`** &mdash; *Function*.



```
  SeisPatch(in::AbstractString,out::AbstractString;<keyword arguments>)
```

Creates overlapping 5d patches from a 5d volume

**Arguments**

  * `in::AbstractString`: input filename (data should have grid information in headers)
  * `out::AbstractString`: prefix for output filenames

**Keyword arguments**

  * `style="sxsygxgy"`: bin style. Options: "mxmyhxhy","mxmyhaz","sxsyhxhy","gxgyhxhy","sxsyhaz","gxgyhaz"
  * `min_isx=0`,`max_isx=0`,`min_isy=0`,`max_isy=0`: grid extreme values for sources
  * `min_igx=0`,`max_igx=0`,`min_igy=0`,`max_igy=0`: grid extreme values for receivers
  * `min_imx=0`,`max_imx=0`,`min_imy=0`,`max_imy=0`: grid extreme values for midpoints
  * `min_ihx=0`,`max_ihx=0`,`min_ihy=0`,`max_ihy=0`: grid extreme values for offsets
  * `min_ih=0`,`max_ih=0`,`min_iaz=0`,`max_iaz=0`: grid extreme values for azimuth and offset
  * `it_WL=9e9`,`it_WO=0` : length and overlapping samples in time patches
  * `ix1_WL=9e9`,`ix1_WO=0`:length and overlapping samples in first space dimension
  * `ix2_WL=9e9`,`ix2_WO=0`,`ix3_WL=9e9`,`ix3_WO=0`,`ix4_WL=9e9`,`ix4_WO=0`

**Output**

`filename,npatch`: String Array with the file name of the data patches, number of patches created

*Credits: A. Stanton, F. Carozzi, 2017*


<a target='_blank' href='https://github.com/fercarozzi/myseismicjulia/tree/c2832f8331d8b4cba573c54c2dd0183c518801d7/src/Utils/SeisPatch.jl#L1-L25' class='documenter-source'>source</a><br>


##SeisUnPatch

<a id='Seismic.SeisUnPatch' href='#Seismic.SeisUnPatch'>#</a>
**`Seismic.SeisUnPatch`** &mdash; *Function*.



```
SeisUnPatch(in,out;<keyword arguments>)
```

Reconstruct a 5D data volume from a set of 5D data patches.

**Arguments**

  * `in::Array{AbstractString,1}`: array containing filename of patches
  * `out::AbstractString`: filename for reconstructed volume

**Keyword arguments**

  * `style="sxsygxgy"`: bin style. Options: "mxmyhxhy","mxmyhaz","sxsyhxhy","gxgyhxhy","sxsyhaz","gxgyhaz"
  * `min_isx=0`,`max_isx=0`,`min_isy=0`,`max_isy=0`: grid extreme values for sources
  * `min_igx=0`,`max_igx=0`,`min_igy=0`,`max_igy=0`: grid extreme values for receivers
  * `min_imx=0`,`max_imx=0`,`min_imy=0`,`max_imy=0`: grid extreme values for midpoints
  * `min_ihx=0`,`max_ihx=0`,`min_ihy=0`,`max_ihy=0`: grid extreme values for offsets
  * `min_ih=0`,`max_ih=0`,`min_iaz=0`,`max_iaz=0`: grid extreme values for azimuth and offset
  * `it_WL=9e9`,`it_WO=0` : length and overlapping samples in time patches
  * `ix1_WL=9e9`,`ix1_WO=0`:length and overlapping samples in first space dimension
  * `ix2_WL=9e9`,`ix2_WO=0`,`ix3_WL=9e9`,`ix3_WO=0`,`ix4_WL=9e9`,`ix4_WO=0`
  * `nt=0`: time samples of reconstructed cube
  * `ang=90`: inline direction measured in degrees CC from East
  * `gamma=1`: vp/vs ratio for PS Asymptotic Conversion Point gathers (use gamma=1 for PP data)
  * `osx=0`,`osy=0`,`ogx=0`,`ogy=0` : origin for source and receiver coordinate system
  * `omx=0`,`omy=0`,`ohx=0`,`ohy=0`: origin for midpoint and offset coordinate system
  * `oaz=0`,`oh=0` : origin for azimuth and offset coordinate system
  * `dsx=1`,`dsy=1`,`dgx=1`,`dgy=1`: source and receiver step-size
  * `dmx=1`,`dmy=1`,`dhx=1`,`dhy=1`: midpoint and offset step-size
  * `dh=1`,`daz=1`: offset and azimuth step-size

**Output**

In file `out`, the 5D reconstructed volume is created.

*Credits: A. Stanton, F Carozzi, 2017*


<a target='_blank' href='https://github.com/fercarozzi/myseismicjulia/tree/c2832f8331d8b4cba573c54c2dd0183c518801d7/src/Utils/SeisUnPatch.jl#L1-L34' class='documenter-source'>source</a><br>


##SeisSort

<a id='Seismic.SeisSort' href='#Seismic.SeisSort'>#</a>
**`Seismic.SeisSort`** &mdash; *Function*.



```
SeisSort(in, out;<keyword arguments>)
```

Sort a seis file using its header words

**Arguments**

  * `in`: input filename >> a text file with information about data extent, data and header file names; a binary file containing data and a binary file containing headers.
  * `out`: output filename

**Keyword arguments**

  * `key=["imx","imy"]`
  * `rev=false` : sort headers in decreasing order
  * `ntrace=1000` : number of traces to read at a time

**Output**

file `out` is created with data sorted.

*Credits: AS, 2015*


<a target='_blank' href='https://github.com/fercarozzi/myseismicjulia/tree/c2832f8331d8b4cba573c54c2dd0183c518801d7/src/Utils/SeisSort.jl#L1-L19' class='documenter-source'>source</a><br>


##SeisWindow

<a id='Seismic.SeisWindow' href='#Seismic.SeisWindow'>#</a>
**`Seismic.SeisWindow`** &mdash; *Function*.



```
SeisWindow(in,out;<keyword arguments>)
```

Window a seis file using header words.

**Arguments**

  * `in::AbstractString`: filename of input
  * `out::AbstractString`: filename of output

**Keyword arguments**

  * `key`
  * `minval`
  * `maxval`

note that windowing along the time axis is achieved by using the key "t".

*Credits: AS, 2015*


<a target='_blank' href='https://github.com/fercarozzi/myseismicjulia/tree/c2832f8331d8b4cba573c54c2dd0183c518801d7/src/Utils/SeisWindow.jl#L1-L18' class='documenter-source'>source</a><br>

