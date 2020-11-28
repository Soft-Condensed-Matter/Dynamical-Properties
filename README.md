Codes to facilitate the computation of time correlations functions in general and 
transport properties in liquids in particular.

The codes are mainly developed in FORTRAN, contributions in other languages are welcome.


Code: MSD.f90

Description: Determines the mean square displacement (MSD.dat) using the coordinates
of N particles contained in the required file MDFps.dat. It is necessary to specify 
the total number of particles and samples within the MDFps.dat file before to compile
the program.


Code: FSelf.f90

Description: Determines the self part of the intermediate scattering function (FSelf.dat)
of N particles contained in the required file MDFps.dat. It is necessary to specify 
the total number of particles and samples within the MDFps.dat file before to compile
the program.


Code: VACF.f90

Description: Determines the velocity autocorrelation function (VACF.dat) of N
particles contained in the required file MDVel.dat. It is necessary to specify 
the total number of particles and samples within the MDVel.dat file before to compile
the program.

