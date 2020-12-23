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


Code: SSCF.f90

Description: Determines the total shear-stress autocorrelation function (TSSCF.dat),
the kynetic (KSSCF.dat), potential (PSSCF.dat) and kynetic-potential (KPSSCF.dat)
contributions of the shear-stress autocorrelation. The program requeries
the instant components given in the files MDTend.dat, MDKtn.dat, MDPtn.dat with the
total, kynetic and potential components, respectively. The number of samples within 
the input files must be provided previously.


Code: Viscocity.f90

Description: Computes the shear-viscocity and its kynetic, potential and 
kynetic-potential contributions using the corresponding autocorrelations functions
given by total shear-stress (TSSCF.dat), kynetic (KSSCF.dat), potential (PSSCF.dat)
and kynetic-potential functions (KPSSCF.dat)


Code: PFACF.f90

Description: Determines the pressure fluctuation autocorrelation function and the 
same correlation of the diagonal tensor pressure component using the file with all 
the components of the tensor pressure given in the file MDTen.dat. The number of 
samples within the input files must be provided previously.


Code: BViscocity.f90

Description: Determines the bulk viscocity from the autocorrelation function of 
the pressure fluctuation by means of the Green-Kubo relation and the longitudal 
viscocity using the fluctuation of the diagonal tensor pressure components. The
program requieres the autoccorelation functions given in the files PFACF.dat and
PFACaa.dat, respectively. 


Code: EFACF.f90

Description: Determines the energy flux (heat current) autocorrelation function 
using its cartesian cartesian components contained in the file MDEfx.dat. 
The number of samples must be previously provided.  
