OUT_STRUCT = confout_mdrun.pdb
GROMPP_FILES = mdout.mdp topol.tpr
MDRUN_FILES = traj.trr ener.edr md.log $(OUT_STRUCT) state.cpt state_prev.cpt traj.xtc traj.tng traj_comp.xtc
SCRIPTS_DIR = ../../scripts

.PHONY: all top clean distclean maintainer-clean

all: confout.gro

conf.pdb:
	$(SCRIPTS_DIR)/gen_ZIF-8_supercell.py 2 2 2 > conf.pdb

zif8_2x2x2_periodic.itp: conf.pdb
	$(SCRIPTS_DIR)/pdb2top_ZIF-8.py conf.pdb > zif8_2x2x2_periodic.itp

confout.gro: $(OUT_STRUCT)
	echo 0 | gmx trjconv -pbc res -f $(OUT_STRUCT) -o confout.gro

$(OUT_STRUCT) traj.trr: topol.tpr
	gmx mdrun -nice 5 -v -c $(OUT_STRUCT)

topol.tpr: grompp.mdp topol.top zif8_2x2x2_periodic.itp
	gmx grompp -pp

top: topol.tpr

#traj_res.xtc: traj.trr
#	echo 0 | gmx trjconv -pbc res -f traj.trr -o traj_res.xtc
#
#angdist.xvg: confout.gro
#	gmx angle -type dihedral -f confout.gro -n angle.ndx -binwidth 1
#
#angle.ndx: zif8_2x2x2_periodic.itp
#	grep "gd_zif8_C2-N-Zn-N" zif8_2x2x2_periodic.itp | awk 'BEGIN { print "[C2-N-Zn-N]" } { print $$1, $$2, $$3, $$4}' | tee $@


clean: distclean

distclean:
	-@rm -f $(GROMPP_FILES) $(MDRUN_FILES)
	-@rm -f \#*\#
	-@rm -f *.xvga energy.xvg processed.top 
	-@rm -f traj_res.xtc angle.ndx angdist.xvg

maintainer-clean: distclean
	-@rm -f conf.pdb zif8_2x2x2_periodic.itp
	-@rm -f confout.gro
