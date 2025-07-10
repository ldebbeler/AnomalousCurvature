Repository: AnomalousCurvature

Features of this project:
* linear radial dispersion, arbitrary exponent 2< alpha < 4 in tangential dispersion
* bare polarization bubble
* real frequencies
* frequency dependence of imaginary part self energy
* frequency and momentum dependence of imaginary part
* momentum dependence of real part self energy via Kramers-Kronig relation

bin-directory: executable binaries
src-directory: source files
include-directory: header files
execute.sh -script: compile and run
data-directory: data in .h5-files and jupyter notebooks for analysis

specalpha.sh is a script for compilation and execution of a program that calculates all quantities for a specific exponent alpha. This parameter is provided as command line argument with the script. The jupyter notebook single_results.ipynb can visualize those results.

alphaRange.sh is a script for compilation and exection of a program that calculates the relevant coefficients of the self energy for a range of parameters alpha. The jupyter notebook coefficients_results.ipynb visualized these results.
