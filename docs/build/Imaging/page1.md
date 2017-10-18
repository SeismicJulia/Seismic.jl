
<a id='Imaging-1'></a>

# Imaging


<a id='PostStackWEM-1'></a>

## PostStackWEM


```
Seismic.PostStackWEM
```


<a id='ShotProfileEWEM-1'></a>

## ShotProfileEWEM

<a id='Seismic.ShotProfileEWEM' href='#Seismic.ShotProfileEWEM'>#</a>
**`Seismic.ShotProfileEWEM`** &mdash; *Function*.



**ShotProfileEWEM**

*Shot Profile Elastic Wave Equation Migration and Demigration of 3D isotropic data.*

**IN**

  * m : vector of filenames of image (comprised of [mpp;mps1;mps2])
  * d : vector of filenames of data (comprised of [ux;uy;uz])
  * adj : flag for adjoint (migration), or forward (demigration) (default=true)
  * damping = 1000. : damping for deconvolution imaging condition
  * vp = "vp" : seis file containing the p-wave velocity (should have same x and z dimensions as the desired image)
  * vs = "vs" : seis file containing the s-wave velocity (should have same x and z dimensions as the desired image)
  * angx = "angx" : seis file containing incidence angles in the x direction for each shot
  * angy = "angy" : seis file containing incidence angles in the y direction for each shot
  * wav = "wav" : seis file containing the source wavelet (in time domain)
  * sz = 0. : source depth (Dev: read this from source wavelet file for variable source depth)
  * gz = 0. : receiver depth (Dev: read this from data file for variable source depth (but then what to do in fwd op?))
  * nangx = 1 : number of angle bins in x direction
  * oangx = 0. : min angle in x direction (angle between source incidence angle and reflector normal in Degrees)
  * dangx = 1. : angle increment in x direction
  * nangy = 1 : number of angle bins in y direction
  * oangy = 0. : min angle in y direction (angle between source incidence angle and reflector normal in Degrees)
  * dangy = 1. : angle increment in y direction
  * fmin = 0. : min frequency to process (Hz)
  * fmax = 80. : max frequency to process (Hz)
  * padt = 1 : pad factor for the time axis
  * padx = 1 : pad factor for the spatial axes
  * verbose = false : flag for error / debugging messages
  * sx = [0.] : array of source X positions (meters)
  * sy = [0.] : array of source Y positions (meters)

**OUT**

*Credits: AS, 2015*


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Imaging/ShotProfileEWEM.jl#L1-L37' class='documenter-source'>source</a><br>


<a id='ShotProfileWEM-1'></a>

## ShotProfileWEM

<a id='Seismic.ShotProfileWEM' href='#Seismic.ShotProfileWEM'>#</a>
**`Seismic.ShotProfileWEM`** &mdash; *Function*.



**ShotProfileWEM**

*Shot Profile Wave Equation Migration and Demigration of 3D isotropic data.*

**IN**

  * m : filename of image
  * d : filename of data
  * adj : flag for adjoint (migration), or forward (demigration) (default=true)
  * pspi : flag for Phase Shift Plus Interpolation (default=true)
  * nref : number of reference velocities to use if pspi is selected (default=5)
  * vel = "vel" : seis file containing the velocity (should have same x and z dimensions as the desired image)
  * angx = "angx" : seis file containing incidence angles in the x direction for each shot
  * angy = "angy" : seis file containing incidence angles in the y direction for each shot
  * wav = "wav" : seis file containing the source wavelet (in time domain)
  * sz = 0. : source depth (Dev: read this from source wavelet file for variable source depth)
  * gz = 0. : receiver depth (Dev: read this from data file for variable source depth (but then what to do in fwd op?))
  * nangx = 1 : number of angle bins in x direction
  * oangx = 0. : min angle in x direction (angle between source incidence angle and reflector normal in Degrees)
  * dangx = 1. : angle increment in x direction
  * nangy = 1 : number of angle bins in y direction
  * oangy = 0. : min angle in y direction (angle between source incidence angle and reflector normal in Degrees)
  * dangy = 1. : angle increment in y direction
  * fmin = 0. : min frequency to process (Hz)
  * fmax = 80. : max frequency to process (Hz)
  * padt = 2 : pad factor for the time axis
  * padx = 2 : pad factor for the spatial axes
  * verbose = false : flag for error / debugging messages
  * sx = [0.] : array of source X positions (meters)
  * sy = [0.] : array of source Y positions (meters)

**OUT**

*Credits: AS, 2015*


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Imaging/ShotProfileWEM.jl#L1-L37' class='documenter-source'>source</a><br>


<a id='ComputeAngles-1'></a>

## ComputeAngles

<a id='Seismic.ComputeAngles' href='#Seismic.ComputeAngles'>#</a>
**`Seismic.ComputeAngles`** &mdash; *Function*.



**ComputeAngles**

*Compute angles for shot gathers. These angles can be used for mapping migrated shots into angle gathers during shot profile migration.*

**IN**

  * angx = "angx" : filename for incidence angles in the x direction for each shot
  * angy = "angy" : filename for incidence angles in the y direction for each shot
  * dip_flag = false : flag to subtract reflector dip from the computed angles to make them with reference to reflector normal
  * vel = "vel" : seis file containing the velocity (should have same x and z dimensions as the desired image)
  * wav = "wav" : seis file containing the source wavelet (in time domain)
  * sz = 0. : source depth (Dev: read this from source wavelet file for variable source depth)
  * nhx = 101 : number of offset bins
  * ohx = 1000. : min offset (surface offset in the data)
  * dhx = 10. : offset increment
  * nhy = 101 : number of offset bins
  * ohy = 1000. : min offset (surface offset in the data)
  * dhy = 10. : offset increment
  * pade_flag = false : flag for Pade Fourier correction
  * fmin = 0. : min frequency to process (Hz)
  * fmax = 80. : max frequency to process (Hz)
  * padt = 2 : pad factor for the time axis
  * padx = 2 : pad factor for the spatial axes
  * verbose = false : flag for error / debugging messages
  * sx = [0.] : array of source X positions (meters)
  * sy = [0.] : array of source Y positions (meters)

**OUT**

*Credits: AS, 2015*


<a target='_blank' href='https://github.com/SeismicJulia/Seismic.jl/blob/42ef65d138b6e379b2d145cd26e18b710f1ae825/src/Imaging/ComputeAngles.jl#L1-L34' class='documenter-source'>source</a><br>

