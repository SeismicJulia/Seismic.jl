# Processing

---

## FKFilterGathers

---

## MatrixMultiply

---

## SeisBandPass

---

## SeisDecimate

---

## SeisFXDecon

**SeisRadonForward**
*Transform a tau-p gather to a time-offset gather using a frequency domain forward Radon operator*
**IN**   
* m: Radon panel,  m[1:nt,1:nh] nt is number of samples and np is number of curvatures or ray parameters
* param={"dt"=>dt,"parab"=>true/false, "p"=>p, "h"=>h, "flow"=>flow, "fhigh"=>fhigh, "href"=>href}
  * dt: sampling interval in secs
  * parab: true --> parabolic transform, false --> linear transform
  * h: offset vector, h[1:nh]
  * if parab=true:  p is a vector of residual moveout at reference offset in seconds (also called curvature)
  * if parab=false: P is a vector of ray parameters  reference (s/m), p[1:np] where np is number of curvatures or ray parameters
  * flow: min frequency in the data (Hz)
  * fhigh: max frequency in the data (Hz)
  * href: reference offset for parabolic Radon Transform
**OUT**  
* d: data synthetized via forward Radon modeling, d[1:nt, 1:nh]

---

## SeisMWNI

---

## SeisNMO

---

## SeisPOCS

---

## SeisPWD

---

## SeisRadon

---

## SeisRadonForward

**SeisRadonForward(m,param=Dict())**

*Transform a tau-p gather to a time-offset gather using a frequency domain forward Radon operator*

**IN**   

* m: Radon panel,  m[1:nt,1:nh] nt is number of samples and np is number of curvatures or ray parameters
* param={"dt"=>dt,"parab"=>true/false, "p"=>p, "h"=>h, "flow"=>flow, "fhigh"=>fhigh, "href"=>href}
  * dt: sampling interval in secs
  * parab: true --> parabolic transform, false --> linear transform
  * h: offset vector, h[1:nh]
  * if parab=true:  p is a vector of residual moveout at reference offset in seconds (also called curvature)
  * if parab=false: P is a vector of ray parameters  reference (s/m), p[1:np] where np is number of curvatures or ray parameters
  * flow: min frequency in the data (Hz)
  * fhigh: max frequency in the data (Hz)
  * href: reference offset for parabolic Radon Transform

**OUT**  

* d: data synthesized via forward Radon modeling, d[1:nt, 1:nh]

*source:*
[Seismic/src/Processing/SeisRadonForward.jl](https://github.com/SeismicJulia/Seismic.jl/tree/b5e44cc4766549fbf044d0040f2c7ef19582b5d2/src/Processing/SeisRadonForward.jl)

---

## SeisRadonInverse

---

## SeisSemblance

---

## SeisStack

---

## SeisWavelets

---

## SmoothGathers

---

## SmoothStructure

---

## fft_op

---
