Imaging
=======

Imaging programs are used to map energy in seismic data to the point of reflection within the subsurface. Their parameters are given by a Dictionary (param = Dict()).

ComputeAngles
^^^^^^^^^^^^^^^
Computation of subsurface source side incidence angles with respect to vertical at all subsurface points for common shot gathers. These angles can be used by ShotProfileWEM to bin migrated shot gathers into angle gathers. A user supplied reflector normal field allows for angles to be corrected to be with reference to reflector normal (see SeisPWD for reflector normal computation).

ShotProfileWEM
^^^^^^^^^^^^^^^
Shot profile wave equation migration or demigration of 3D acoustic data using an isotropic velocity model. Migrated shot gathers can be binned into 3D angle gathers (angx, angy) by bilinear interpolation into angle bins, where angles for each shot gather are user supplied (see ComputeAngles for angle compuation).

ShotProfileLSWEM
^^^^^^^^^^^^^^^
Least squares migration using ShotProfileWEM as a linear operator.

