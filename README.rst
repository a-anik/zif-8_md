======================
ZIF-8 GROMACS MD files
======================

Files for molecular dynamics simulation of ZIF-8 metal-organic framework with embeded guest molecules (TEMPO nitroxide radical, CO2, N2) [1]_.

Prerequisites
-------------
* GROMACS 2016.3 molecular dynamics package
* GNU make
* coreutils, sed
* python 2.7, numpy, scipy (for some analyis scripts)
* MDAnalysis python library (for ZIF-8 topology generator)


Directory contents
------------------
* scripts/ - common analyis scripts, makefiles, simulation parameters;
* top_all/ - molecular topologies;
* ZIF-8_TEMPO_03co2/ - simulation template dir for 2x2x2 ZIF-8 supercell with 1 TEMPO and 24 CO2 (3 per u.c.) molecules;
* ZIF-8_TEMPO_Ngas_template/  - general template dir for ZIF-8 + TEMPO + GAS simulations.

Citation
--------
.. [1] A.M. Sheveleva, A.V. Anikeenko, A.S. Poryvaev, D.L. Kuzmina, I.K. Shundrina, D.I. Kolokolov, A.G. Stepanov, M.V. Fedin, Probing Gas Adsorption in Metal–Organic Framework ZIF-8 by EPR of Embedded Nitroxides, J. Phys. Chem. C. 121 (2017) 19880–19886. `doi:10.1021/acs.jpcc.7b06884 <http://dx.doi.org/10.1021/acs.jpcc.7b06884>`_.
