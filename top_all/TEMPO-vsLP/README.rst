================================
TEMPO parametrization from [1]_.
================================

Changes 
-------

- AMBER topology was converted to GROMACS with ACPYPE [2]_
- dihedral angles of multiple-dihedral-type were used instead of Ryckaert-Bellemans (gmx45 option)
- topology was divided into separate files (atomtypes and molecule .itp)
- added suffix to TEMPO's atom types
- lone pairs LP1, LP2 were converted to massless virtual sites ([virtual_sites3] type in GROMACS).


.. [1] \E. Stendardo, A. Pedone, P. Cimino, M. C. Menziani, O. Crescenzi, and V. Barone, “Extension of the AMBER force-field for the study of large nitroxides in condensed phases: an ab initio parameterization,” Phys. Chem. Chem. Phys., vol. 12, no. 37, pp. 11697–11709, Sep. 2010.

.. [2] \A. W. Sousa da Silva and W. F. Vranken, “ACPYPE - AnteChamber PYthon Parser interfacE,” BMC Res Notes, vol. 5, p. 367, Jul. 2012.

