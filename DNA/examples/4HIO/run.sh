#! /usr/bin/bash
mkdir Output
echo "Using nucproSBM.py, a part to gokit for generating CG GROMACS input files for nucleotide protein system"
echo "Using 4HIO.pdb as the test input."
echo "Help output have been generated in USAGE.txt"
python ../../nucproSBM.py --help >> USAGE.txt
echo "======================================================="
echo "Case I: 4HIO.pdb as input pdb. Using charge on CB (centre at the co-ordinate of farthest atom) with Monovalent ionic conc to be 0.01M. CA and CB radius set to 1.9 and 1.5A respectively. For DNA beads, using default radius (refer to USAGE.txt). Using Debye-Huckel electrostatic potential"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain"

echo "Ready to proceed [yes=1/no=0]? (0)"
read input
if [ $input -eq 1 ]
then
python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain 

mkdir Output/CaseI
rm -r Output/CaseI/*
mv * Output/CaseI/
cp Output/CaseI/4HIO.pdb .
cp Output/CaseI/run.sh .
mv Output/CaseI/polyT.pdb .
mv Output/CaseI/4HIO_DNA.pdb .
fi

echo "======================================================="

echo "Case II: CB at Sidechain COM"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro  --skip_glycine  --CBcharge --dielec 70 --iconc 0.01 --debye --all_chain"
echo "Ready to proceed [yes=1/no=0]?(0)"
read input
if [ $input -eq 1 ]
then
python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro   --skip_glycine  --CBcharge --dielec 70 --iconc 0.01 --debye --all_chain
mkdir Output/CaseII
rm -r Output/CaseII/*
mv * Output/CaseII/
cp Output/CaseII/4HIO.pdb .
cp Output/CaseII/run.sh .
mv Output/CaseII/polyT.pdb .
mv Output/CaseII/4HIO_DNA.pdb .
fi
python ../../nucproSBM.py --help >> USAGE.txt

echo "======================================================="
echo "Case III: 4HIO.pdb as input pdb. Using different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface"

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

python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface #--custom_nuc polyT.pdb

mkdir Output/CaseIII
rm -r Output/CaseIII/*
mv * Output/CaseIII/
cp Output/CaseIII/4HIO.pdb .
cp Output/CaseIII/run.sh .
mv Output/CaseIII/polyT.pdb .
mv Output/CaseIII/4HIO_DNA.pdb .
fi

echo "======================================================="
echo "Case IV: 4HIO.pdb as input pdb. Using different DNA structure along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. The file is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "The DNA sequence used is PolyT of lenght 9 bases (same as the native bound DNA strand)"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb"
cp ../polyT.pdb .
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

python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb

mkdir Output/CaseIV
rm -r Output/CaseIV/*
mv * Output/CaseIV/
cp Output/CaseIV/4HIO.pdb .
cp Output/CaseIV/polyT.pdb .
cp Output/CaseIV/run.sh .
mv Output/CaseIV/polyT.pdb .
mv Output/CaseIV/4HIO_DNA.pdb .
fi

echo "======================================================="
echo "Case V: 4HIO.pdb as input pdb. Using native DNA sequecne structure as custom input (no notive contacts) along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. The file is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 4HIO_DNA.pdb"
cp ../4HIO_DNA.pdb .
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

python ../../nucproSBM.py --aa_pd 4HIO.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 4HIO_DNA.pdb

mkdir Output/CaseV
rm -r Output/CaseV/*
mv * Output/CaseV/
cp Output/CaseV/4HIO.pdb .
cp Output/CaseV/4HIO_DNA.pdb .
cp Output/CaseV/run.sh .
mv Output/CaseV/polyT.pdb .
mv Output/CaseV/4HIO_DNA.pdb .
fi
