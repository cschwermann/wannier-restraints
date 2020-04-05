# Wannier-Restraints
The scripts / programs in this repository are useful for the analysis of results produced by the modified CPMD (https://www.cpmd.org) code presented in
C. Schwermann and N. L. Doltsinis. “Exciton transfer free energy from Car–Parrinello molecular dynamics”. Phys. Chem. Chem. Phys. (2020). Advance Article (https://doi.org/10.1039/C9CP06419B).

## cpwan
cpwan is a short script which renames the "WC_POT.xyz" trajectories from several subfolders, arranged by their reaction coordinate, to WC_POT_1.xyz, WC_POT_2.xyz, etc. The script copies these files either to the current directory, or if an argument is given, to the specified directory. The renamed files can then be used as input for the umbrella.f90 analysis program.

Example: three calculations with Wannier-Restraints for reaction coordiantes 0.0, 0.5 and 1.0 were performed and the results are present in the subfolders 0.0/, 0.5/ and 1.0/. The command 
    `cpwan run1`
creates the folder run1/, which then contains the three files WC_POT_1.xyz, WC_POT_2.xyz and WC_POT_3.xyz, corresponding to the results from the folders 0.0/, 0.5/ and 1.0/.

cpwan might only really be useful if this folder structure and naming scheme is used...

## umbrella.f90
umbrella.f90 is an analysis tool which calculates the free energy of exciton transfer from the WC_POT.xyz trajectories from CPMD calculations using Wannier-Restraints.

The program expects <nwind> files to be present in the directory, named &lt;filnam>_&lt;i>.xyz, with &lt;i> in [1,&lt;nwind>].

Usage:
    `umbrella <nwind> <temp> <xi_start> <xi_end> <orb. pos.> <filename> [#integration points]`
Arguments in <> are necessary, arguments in [ ] are optional.
* &lt;nwind> is the number of sampling windows,
* &lt;temp> is the simulation temperature,
* &lt;xi_start> and &lt;xi_end> are start and end point of the reaction coordinate,
* &lt;orb. pos.> is the position of the orbital of interest in the xyz files, 
* &lt;filename> is the filename (see above),
* [#integration points] is the number of integration points for the free energy, default is 1000.

The program prints a list containing the number of the window, the average reaction coordinate for that window, the standard deviation of the reaction coordinate, the number of steps in the respective trajectory, the force constant of the restraint and the reference coordinate. Additionally, a file "Fenergy.dat" is produced, containing the calculated free energy profile. The first column in that file is the reaction coordiante (dimensionless), the second column is the free energy (in electronvolt) and the third column contains the average restraint force (in electronvolt per reaction coordinate).

umbrella.f90 can be compiled with any Fortran  compiler, e.g.
    `gfortran -o umbrella umbrella.f90`
