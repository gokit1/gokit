gmx grompp -f md.mdp  -c gromacs.gro  -p gromacs.top  -po mdout.mdp -o run.tpr -maxwarn 1 > a 2> b



gmx mdrun -x traj.xtc -e ener.edr -o traj.trr -s run.tpr -g md.log -table table_file.xvg > run.out 2> run.err
