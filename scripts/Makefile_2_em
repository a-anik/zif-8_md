OUT_STRUCT = confout.gro
GROMPP_FILES = mdout.mdp topol.tpr
MDRUN_FILES = traj.trr ener.edr md.log $(OUT_STRUCT) state.cpt state_prev.cpt traj.xtc traj_comp.xtc

.PHONY: clean all

all: run

run: top
	gmx mdrun -nice 5 -v -c $(OUT_STRUCT)

topol.tpr: grompp.mdp topol.top
	gmx grompp

top: topol.tpr


clean:
	-@rm -f $(GROMPP_FILES) $(MDRUN_FILES)
	-@rm -f \#*\#
	-@rm -f energy.xvg
