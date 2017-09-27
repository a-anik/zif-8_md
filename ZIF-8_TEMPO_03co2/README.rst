===========
Description
===========
Sample directory for the molecular dynamics simulation of ZIF-8 2x2x2 supercell + 1 TEMPO molecule + 24 CO2 molecules.

Usage
-----
- type "make all" in the toplevel directory

Steps
-----
1. Insertion of guest molecules into ZIF-8
2. Energy minimization
3. Equilibration : 1 ns, 300 K, NVT
4. Production run : 100 ns, 300 K, NVT (only ZIF-8 atoms are coupled to the thermal bath)

Additional analysis steps
-------------------------
5. Calculation of Rotational autocorrelation funcions (for TEMPO molecule)
6. Calculation of Spatial distribution functions (CO2, TEMPO atoms)
7. Calculation of Distance distribution functions, counting number of CO2 molecules in central cage.

Run produces around 12Gb of data (requires ~50Gb in the process), takes ~9h on i7-4790K CPU with GeForce GTX 980 GPU accelerator.
