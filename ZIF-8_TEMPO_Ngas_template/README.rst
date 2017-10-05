===========
Description
===========
Template directory for the molecular dynamics simulation of 2x2x2 supercell of ZIF-8 MOF + 1 TEMPO molecule + N gas molecules.

Usage
-----

- edit gas.mk file to choose between CO2/N2 and set the number of gas molecules
- type ``make all`` in the toplevel directory

Simulation steps
----------------

1. Random insertion of guest molecules into ZIF-8
2. Energy minimization
3. Equilibration: 1 ns, NVT, 300 K
4. Production run: 100 ns, 300 K, NVT (only ZIF-8 atoms are coupled to the thermal bath)


Run produces around 12Gb of data (requires ~50Gb in process), takes ~9h on i7-4790K CPU @ 4.00GHz with GeForce GTX 980 acceleration.
