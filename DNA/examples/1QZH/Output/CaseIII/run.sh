#! /usr/bin/bash
mkdir Output
echo "Using nucproSBM.py, a part to gokit for generating CG GROMACS input files for nucleotide protein system"
echo "Using 1QZH.pdb as the test input."
echo "Help output have been generated in USAGE.txt"
python ../../nucproSBM.py --help >> USAGE.txt
echo "======================================================="
echo "Case I: 1QZH.pdb as input pdb. Using charge on CB (centre at the co-ordinate of farthest atom) with Monovalent ionic conc to be 0.01M. CA and CB radius set to 1.9 and 1.5A respectively. For DNA beads, using default radius (refer to USAGE.txt). Using Debye-Huckel electrostatic potential"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 1QZH.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye"

echo "Ready to proceed [yes=1/no=0]? (0)"
read input
if [ $input -eq 1 ]
then
python ../../nucproSBM.py --aa_pd 1QZH.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye

mkdir Output/CaseI
rm -r Output/CaseI/*
mv * Output/CaseI/
cp Output/CaseI/1QZH.pdb .
cp Output/CaseI/run.sh .
fi
echo "======================================================="

echo "Case II: Default radius of all Beads"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 1QZH.pdb --grotop gromacs.top -pdbgro gromacs.gro  --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --debye"
echo "Ready to proceed [yes=1/no=0]?(0)"
read input
if [ $input -eq 1 ]
then
python ../../nucproSBM.py --aa_pd 1QZH.pdb --grotop gromacs.top -pdbgro gromacs.gro   --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --debye
mkdir Output/CaseII
rm -r Output/CaseII/*
mv * Output/CaseII/
cp Output/CaseII/1QZH.pdb .
cp Output/CaseII/run.sh .
fi
python ../../nucproSBM.py --help >> USAGE.txt

echo "======================================================="
echo "Case III: 1QZH.pdb as input pdb. Using different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 1QZH.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --interface"

echo "Ready to proceed [yes=1/no=0]? (0)"
read input
if [ $input -eq 1 ]
then
bash ../../generate_NucproInterface_inputfile.sh
echo "

 2.0
 2.0
 2.0
 2.0
 
 2.0
 2.0



" > myparams.txt
grep -v "^;" nucpro_interface.input > temp1.txt
paste temp1.txt myparams.txt > temp2.txt
mv temp2.txt nucpro_interface.input
rm myparams.txt temp1.txt temp2.txt
python ../../nucproSBM.py --aa_pd 1QZH.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --interface

mkdir Output/CaseIII
rm -r Output/CaseIII/*
mv * Output/CaseIII/
cp Output/CaseIII/1QZH.pdb .
cp Output/CaseIII/run.sh .
fi
echo "======================================================="

echo "Case IV: NOT removing unique chains (Keep all chains)"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 1QZH.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chains" 
echo "Warning: This will take time"
echo "Ready to proceed [yes=1/no=0]? (0)"
read input
if [ $input -eq 1 ]
then
#echo "Start_time" >> ../1QZH_all_chain_runtime.log
#date >> ../1QZH_all_chain_runtime.log
python ../../nucproSBM.py --aa_pd 1QZH.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chains
#echo "End_time" >> ../1QZH_all_chain_runtime.log
#date >> ../1QZH_all_chain_runtime.log
mkdir Output/CaseIV
rm -r Output/CaseIV/*
mv * Output/CaseIV/
cp Output/CaseIV/1QZH.pdb .
cp Output/CaseIV/run.sh .
fi

