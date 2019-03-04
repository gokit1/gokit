
grompp -f md.mdp -c gromacs.gro -p gromacs.top -po mdout.out -o run.tpr -maxwarn 1 > a 2>b
mdrun -x traj.xtc -e ener.edr -o traj.trr -s run.tpr -g md.log -table table_file.xvg -tablep table_file.xvg > run.out 2> run.err
