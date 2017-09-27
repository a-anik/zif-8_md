OUT_STRUCT = confout.pdb
GROMPP_FILES = mdout.mdp topol.tpr
MDRUN_FILES = traj.trr ener.edr md.log $(OUT_STRUCT) state.cpt state_prev.cpt traj.xtc traj_comp.xtc

.PHONY: clean all top distclean maintainer-clean


all: TEMPO-vsLP_centered.pdb

TEMPO-vsLP_centered.pdb: $(OUT_STRUCT)
	gmx editconf -center 0 -f $(OUT_STRUCT) -o $@

$(OUT_STRUCT): topol.tpr
	gmx mdrun -nice 5 -v -c $(OUT_STRUCT)

topol.tpr: grompp.mdp topol.top
	gmx grompp

top: topol.tpr

clean: distclean

distclean:
	-@rm -f $(GROMPP_FILES) $(MDRUN_FILES)
	-@rm -f \#*\#
	-@rm -f energy.xvg res_trajout.xtc rmsd.xvg

maintainer-clean: distclean
	-@rm -f TEMPO-vsLP_centered.pdb
