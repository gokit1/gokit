#! /usr/bin/bash
mkdir Output
echo "Using nucproSBM.py, a part to gokit for generating CG GROMACS input files for nucleotide protein system"
echo "Using 1QYS.pdb as the test input."
echo "1QZH (Top7) is an example of PROTEIN ONLY case. In such cases all the paramters corrosponding to DNA/RNA will be ignored. The output will be same as the protein only code."
echo "Although nucproSBM.py can handle protein only cases, for more options its better to use the protein only code"
echo "Help output have been generated in USAGE.txt"
python ../../nucproSBM.py --help >> USAGE.txt
echo "======================================================="
echo "Case I: 1qys_H.pdb as input pdb. Using charge on CB (centre at the co-ordinate of farthest atom) with Monovalent ionic conc to be 0.01M. CA and CB radius set to 1.9 and 1.5A respectively. For DNA beads, using default radius (refer to USAGE.txt). Using Debye-Huckel electrostatic potential"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 1qys_H.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye"

echo "Ready to proceed [yes=1/no=0]? (0)"
read input
if [ $input -eq 1 ]
then
python ../../nucproSBM.py --aa_pd 1qys_H.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye

mkdir Output/CaseI
rm -r Output/CaseI/*
mv * Output/CaseI/
cp Output/CaseI/1qys_H.pdb .
cp Output/CaseI/run.sh .
fi
echo "======================================================="
python ../../nucproSBM.py --help >> USAGE.txt