# Reading and Writing Seismic data from disk

Seismic.jl has the ability to read and write a number of standard seismic data formats and
uses its own binary representation of seismic data on disk that is based on splitting the 
file into data (filename.seisd) and headers (filename.seish).

---

## The .seisd/.seish format

The Seismic.jl package uses a binary representation of seismic data files on disk.
A seismic file converted to the Seismic.jl format is broken into two parts:
`foo.seisd` containing seismic data as IEEE 32-bit floats, and `foo.seish` 
containing 31 32bit headers for each trace.

To see the header that are contained in the type `Header`, you can type:

```no-highlight
julia> names(Header)
31-element Array{Symbol,1}:
 :tracenum
 :o1      
 :n1      
 :d1      
 :sx      
 :sy      
 :gx      
 :gy      
 :mx      
 :my      
 :hx      
 :hy      
 :h       
 :az      
 :ang     
 :isx     
 :isy     
 :igx     
 :igy     
 :imx     
 :imy     
 :ihx     
 :ihy     
 :ih      
 :iaz     
 :iang    
 :selev   
 :gelev   
 :sstat   
 :gstat   
 :trid
```

## Reading SEGY, SeismicUnix, and RSF (Madagascar) files


A simple example of file conversion looks like this:

```no-highlight
SegyToSeis("foo.su","foo",{"format"=>su})
```

which would create `foo.seisd` and `foo.seish`.

## Writing SeismicUnix files

Converting from a .seisd/.seish file to SeismicUnix format

```no-highlight
SeisToSegy("foo","foo.su",{"format"=>su})
```

which would create `foo.su`.

## Reading data and headers

A simple example of reading a .seisd/.seish file looks like:

```no-highlight
d,h = SeisRead("foo")
```

Here we read the files foo.seisd and put it into the variable d. The type of d is Array{Float32,2)
with dimensions number of time samples x number of traces. We also read the file foo.seish and put
it into the variable h. The type of h is Array{Header,1} (a vector with elements of type Header). 

## Reading headers only

If we only want to read the headers we can use the function `SeisReadHeaders`

```no-highlight
h = SeisReadHeaders("foo")
```

Here we read the file foo.seish and put it into the variable h. 

## Inspecting the Headers of a .seish file

To get some statistics for a .seish file you can use the command `SeisHeaderInfo`. For
example, the headers for the open source Teapot Dome 3D survey are found to be:

```no-highlight
julia> SeisHeaderInfo("npr3_gathers_nmo")
Displaying information for npr3_gathers_nmo.seish (723991 traces):
       Key          Minimum          Maximum             Mean
=============================================================
  tracenum            1.000       723991.000       361996.000
        o1            0.000            0.000            0.000
        n1         2049.000         2049.000         2049.000
        d1            0.002            0.002            0.002
        sx       788339.938       809248.625       800467.813
        sy       939328.375       976916.125       959490.438
        gx       788173.000       809316.500       800253.938
        gy       939209.813       976867.938       959768.188
        mx       788256.500       809282.563       800360.938
        my       939269.125       976892.000       959629.313
        hx       -18298.375        18035.438         -213.839
        hy       -10737.063        15758.063          277.830
         h          100.222        19276.408         6224.425
        az            0.000          359.999          186.337
       ang            0.000            0.000            0.000
       isx       788340.000       809249.000       800467.813
       isy       939328.000       976916.000       959490.438
       igx       788173.000       809317.000       800254.125
       igy       939210.000       976868.000       959768.188
       imx           -6.000          185.000          103.857
       imy            4.000          346.000          188.951
       ihx       -18298.000        18035.000         -213.840
       ihy       -10737.000        15758.000          277.830
        ih            0.000           46.000           14.820
       iaz            0.000            8.000            4.140
      iang            0.000            0.000            0.000
     selev         4934.400         5405.800         5158.012
     gelev         4931.300         5656.800         5160.850
     sstat            0.000            0.000            0.000
     gstat            0.000            0.000            0.000
      trid            0.000            2.000            1.003

```

## Extracting a Header word

To extract a single header word from a vector of headers you can use the internal command `Seismic.ExtractHeader`.
For example, to plot source and receiver coordinates for the Teapot Dome dataset we can type:

```no-highlight
julia> Pkg.add("ASCIIPlots"); using ASCIIPlots;
julia> h = SeisReadHeaders("npr3_gathers_nmo");
julia> sx = Seismic.ExtractHeader(h,"sx");
julia> sy = Seismic.ExtractHeader(h,"sy");
julia> gx = Seismic.ExtractHeader(h,"gx");
julia> gy = Seismic.ExtractHeader(h,"gy");
julia> scatterplot(sx,sy,sym='X')

	-------------------------------------------------------------
	|                      X                                     | 976916.13
	|          XXXXXX XXXXXXXXXXX XXX                            |
	|      XXXXXXXXXXXXXXXXXX XXXXXXXXXX                         |
	|       XXXXX XXXXX XXXXXXXXXXXXXXXXX XXX                    |
	|   XXXXXXXXXX XXXXX XXXXXXXXXXX XXXXX XX                    |
	|XX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXX                 |
	|      XXXXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXX  X             |
	|           XXXXXXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXX            |
	|               X  XXXX XXXXXXXXXXX XXXXX XXXXXXXXXXX        |
	|                   XXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXX XX     |
	|                      XXXXXXXXXXXXXXXXXXX XXXX X XXXXXX     |
	|                          XXXX XXXXXXXXXXXXXX XXXXXXXXX     |
	|                        XXXXXXXXXXXXXXXXXXX XXXXXXXXXXX X   |
	|                          XXXXXXXXXXXX XXXXX XXXXXXXXXXXXXX |
	|                           X  XXXXXXXXXXXXXXXX XXXXXXXXXXXX |
	|                              XXX XXXXXXXXXXXXXXXXXXXXXXXXXX|
	|                                 X  XXXX XXXXX XXXXXXXXXXX  |
	|                                      XXXXXXXXX XXXXXXXXXX  |
	|                                         X   XXX XXXXXX     |
	|                                             XXXX X         | 939328.38
	-------------------------------------------------------------
	788339.94                                                    809248.63

julia> scatterplot(gx,gy,sym='^')

	-------------------------------------------------------------
	|                            ^                               | 976867.94
	|          ^^^^^^^^^^^^^^^^^^^^^^^                           |
	|       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                        |
	|   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                    |
	|   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                |
	|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                |
	|   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^             |
	|           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^             |
	|               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         |
	|                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     |
	|                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     |
	|                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     |
	|                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ |
	|                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ |
	|                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ |
	|                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|
	|                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ |
	|                                  ^^^^^^^^^^^^^^^^^^^^^^^^^ |
	|                                         ^^^^^^^^^^^^^^^^^^ |
	|                                         ^^^^^^^^^^^^^^^^^^ | 939209.81
	-------------------------------------------------------------
	788173.00                                                    809316.50

```






