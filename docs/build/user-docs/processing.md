
<a id='Processing-1'></a>

# Processing


---


<a id='FKFilterGathers-1'></a>

## FKFilterGathers


---


<a id='MatrixMultiply-1'></a>

## MatrixMultiply


---


<a id='SeisBandPass-1'></a>

## SeisBandPass


---


<a id='SeisDecimate-1'></a>

## SeisDecimate


---


<a id='SeisFXDecon-1'></a>

## SeisFXDecon


**SeisRadonForward** *Transform a tau-p gather to a time-offset gather using a frequency domain forward Radon operator* **IN**   


  * m: Radon panel,  m[1:nt,1:nh] nt is number of samples and np is number of curvatures or ray parameters
  * param={"dt"=>dt,"parab"=>true/false, "p"=>p, "h"=>h, "flow"=>flow, "fhigh"=>fhigh, "href"=>href}

      * dt: sampling interval in secs
      * parab: true –> parabolic transform, false –> linear transform
      * h: offset vector, h[1:nh]
      * if parab=true:  p is a vector of residual moveout at reference offset in seconds (also called curvature)
      * if parab=false: P is a vector of ray parameters  reference (s/m), p[1:np] where np is number of curvatures or ray parameters
      * flow: min frequency in the data (Hz)
      * fhigh: max frequency in the data (Hz)
      * href: reference offset for parabolic Radon Transform


**OUT**  


  * d: data synthetized via forward Radon modeling, d[1:nt, 1:nh]


---


<a id='SeisMWNI-1'></a>

## SeisMWNI


---


<a id='SeisNMO-1'></a>

## SeisNMO


---


<a id='SeisPOCS-1'></a>

## SeisPOCS


---


<a id='SeisPWD-1'></a>

## SeisPWD


---


<a id='SeisRadon-1'></a>

## SeisRadon


---


<a id='SeisRadonForward-1'></a>

## SeisRadonForward


**SeisRadonForward(m,param=Dict())**


*Transform a tau-p gather to a time-offset gather using a frequency domain forward Radon operator*


**IN**   


  * m: Radon panel,  m[1:nt,1:nh] nt is number of samples and np is number of curvatures or ray parameters
  * param={"dt"=>dt,"parab"=>true/false, "p"=>p, "h"=>h, "flow"=>flow, "fhigh"=>fhigh, "href"=>href}

      * dt: sampling interval in secs
      * parab: true –> parabolic transform, false –> linear transform
      * h: offset vector, h[1:nh]
      * if parab=true:  p is a vector of residual moveout at reference offset in seconds (also called curvature)
      * if parab=false: P is a vector of ray parameters  reference (s/m), p[1:np] where np is number of curvatures or ray parameters
      * flow: min frequency in the data (Hz)
      * fhigh: max frequency in the data (Hz)
      * href: reference offset for parabolic Radon Transform


**OUT**  


  * d: data synthesized via forward Radon modeling, d[1:nt, 1:nh]


*source:* [Seismic/src/Processing/SeisRadonForward.jl](https://github.com/SeismicJulia/Seismic.jl/tree/b5e44cc4766549fbf044d0040f2c7ef19582b5d2/src/Processing/SeisRadonForward.jl)


---


<a id='SeisRadonInverse-1'></a>

## SeisRadonInverse


---


<a id='SeisSemblance-1'></a>

## SeisSemblance


---


<a id='SeisStack-1'></a>

## SeisStack


---


<a id='SeisWavelets-1'></a>

## SeisWavelets


---


<a id='SmoothGathers-1'></a>

## SmoothGathers


---


<a id='SmoothStructure-1'></a>

## SmoothStructure


---


<a id='fft_op-1'></a>

## fft_op


---

