mkdir Output
mkdir Output/Charged_Aromatic
mkdir Output/Charged
mkdir Output/Aromatic
mkdir Output/Blank_control
echo "Creating files for interace Charge + Aromatic"
echo "======================================================="
echo "nat_DNA_nat_contacts: 2C62.pdb as input pdb. Using different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface"

echo "Ready to proceed [yes=1/no=0]? (0)"
#read input
input=1
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface #--custom_nuc polyT.pdb

mkdir Output/Charged_Aromatic/nat_DNA_nat_contacts
rm -r Output/Charged_Aromatic/nat_DNA_nat_contacts/*
mv * Output/Charged_Aromatic/nat_DNA_nat_contacts/
cp Output/Charged_Aromatic/nat_DNA_nat_contacts/2C62.pdb .
cp Output/Charged_Aromatic/nat_DNA_nat_contacts/run.sh .
mv Output/Charged_Aromatic/nat_DNA_nat_contacts/polyT.pdb .
mv Output/Charged_Aromatic/nat_DNA_nat_contacts/2C62_DNA.pdb .
fi

echo "======================================================="
echo "nonat_DNA_nonat_contacts: 2C62.pdb as input pdb. Using different DNA structure along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "The DNA sequence used is PolyT of lenght 11bases (same as the native bound DNA strand)"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb"
cp ../polyT.pdb .
echo "Ready to proceed [yes=1/no=0]? (0)"
#read input
input=0
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb

mkdir Output/Charged_Aromatic/nonat_DNA_nonat_contacts
rm -r Output/Charged_Aromatic/nonat_DNA_nonat_contacts/*
mv * Output/Charged_Aromatic/nonat_DNA_nonat_contacts/
cp Output/Charged_Aromatic/nonat_DNA_nonat_contacts/2C62.pdb .
cp Output/Charged_Aromatic/nonat_DNA_nonat_contacts/polyT.pdb .
cp Output/Charged_Aromatic/nonat_DNA_nonat_contacts/run.sh .
mv Output/Charged_Aromatic/nonat_DNA_nonat_contacts/polyT.pdb .
mv Output/Charged_Aromatic/nonat_DNA_nonat_contacts/2C62_DNA.pdb .
fi

echo "======================================================="
echo "nat_DNA_nonat_contacts: 2C62.pdb as input pdb. Using native DNA sequecne structure as custom input (no notive contacts) along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 2C62_DNA.pdb"
cp ../2C62_DNA.pdb .
echo "Ready to proceed [yes=1/no=0]? (0)"
#read input
input=0
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 2C62_DNA.pdb

mkdir Output/Charged_Aromatic/nat_DNA_nonat_contacts
rm -r Output/Charged_Aromatic/nat_DNA_nonat_contacts/*
mv * Output/Charged_Aromatic/nat_DNA_nonat_contacts/
cp Output/Charged_Aromatic/nat_DNA_nonat_contacts/2C62.pdb .
cp Output/Charged_Aromatic/nat_DNA_nonat_contacts/2C62_DNA.pdb .
cp Output/Charged_Aromatic/nat_DNA_nonat_contacts/run.sh .
mv Output/Charged_Aromatic/nat_DNA_nonat_contacts/polyT.pdb .
mv Output/Charged_Aromatic/nat_DNA_nonat_contacts/2C62_DNA.pdb .
fi
echo "Charged interface contacts only"
echo "======================================================="
echo "nat_DNA_nat_contacts: 2C62.pdb as input pdb. Using different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface"

echo "Ready to proceed [yes=1/no=0]? (0)"
#read input
input=0
if [ $input -eq 1 ]
then
bash ../../generate_NucproInterface_inputfile.sh
echo "
False
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface #--custom_nuc polyT.pdb

mkdir Output/Charged/nat_DNA_nat_contacts
rm -r Output/Charged/nat_DNA_nat_contacts/*
mv * Output/Charged/nat_DNA_nat_contacts/
cp Output/Charged/nat_DNA_nat_contacts/2C62.pdb .
cp Output/Charged/nat_DNA_nat_contacts/run.sh .
mv Output/Charged/nat_DNA_nat_contacts/polyT.pdb .
mv Output/Charged/nat_DNA_nat_contacts/2C62_DNA.pdb .
fi

echo "======================================================="
echo "nonat_DNA_nonat_contacts: 2C62.pdb as input pdb. Using different DNA structure along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "The DNA sequence used is PolyT of lenght 11bases (same as the native bound DNA strand)"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb"
cp ../polyT.pdb .
echo "Ready to proceed [yes=1/no=0]? (0)"
#read input
input=0
if [ $input -eq 1 ]
then
bash ../../generate_NucproInterface_inputfile.sh
echo "
False
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb

mkdir Output/Charged/nonat_DNA_nonat_contacts
rm -r Output/Charged/nonat_DNA_nonat_contacts/*
mv * Output/Charged/nonat_DNA_nonat_contacts/
cp Output/Charged/nonat_DNA_nonat_contacts/2C62.pdb .
cp Output/Charged/nonat_DNA_nonat_contacts/polyT.pdb .
cp Output/Charged/nonat_DNA_nonat_contacts/run.sh .
mv Output/Charged/nonat_DNA_nonat_contacts/polyT.pdb .
mv Output/Charged/nonat_DNA_nonat_contacts/2C62_DNA.pdb .
fi

echo "======================================================="
echo "nat_DNA_nonat_contacts: 2C62.pdb as input pdb. Using native DNA sequecne structure as custom input (no notive contacts) along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 2C62_DNA.pdb"
cp ../2C62_DNA.pdb .
echo "Ready to proceed [yes=1/no=0]? (0)"
#read input
input=0
if [ $input -eq 1 ]
then
bash ../../generate_NucproInterface_inputfile.sh
echo "
False
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 2C62_DNA.pdb

mkdir Output/Charged/nat_DNA_nonat_contacts
rm -r Output/Charged/nat_DNA_nonat_contacts/*
mv * Output/Charged/nat_DNA_nonat_contacts/
cp Output/Charged/nat_DNA_nonat_contacts/2C62.pdb .
cp Output/Charged/nat_DNA_nonat_contacts/2C62_DNA.pdb .
cp Output/Charged/nat_DNA_nonat_contacts/run.sh .
mv Output/Charged/nat_DNA_nonat_contacts/polyT.pdb .
mv Output/Charged/nat_DNA_nonat_contacts/2C62_DNA.pdb .
fi
echo "Aromatic cross terms only"
echo "======================================================="
echo "nat_DNA_nat_contacts: 2C62.pdb as input pdb. Using different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface"

echo "Ready to proceed [yes=1/no=0]? (0)"
#read input
input=0
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface #--custom_nuc polyT.pdb

mkdir Output/Aromatic/nat_DNA_nat_contacts
rm -r Output/Aromatic/nat_DNA_nat_contacts/*
mv * Output/Aromatic/nat_DNA_nat_contacts/
cp Output/Aromatic/nat_DNA_nat_contacts/2C62.pdb .
cp Output/Aromatic/nat_DNA_nat_contacts/run.sh .
mv Output/Aromatic/nat_DNA_nat_contacts/polyT.pdb .
mv Output/Aromatic/nat_DNA_nat_contacts/2C62_DNA.pdb .
fi

echo "======================================================="
echo "nonat_DNA_nonat_contacts: 2C62.pdb as input pdb. Using different DNA structure along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "The DNA sequence used is PolyT of lenght 11bases (same as the native bound DNA strand)"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb"
cp ../polyT.pdb .
echo "Ready to proceed [yes=1/no=0]? (0)"
#read input
input=0
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb

mkdir Output/Aromatic/nonat_DNA_nonat_contacts
rm -r Output/Aromatic/nonat_DNA_nonat_contacts/*
mv * Output/Aromatic/nonat_DNA_nonat_contacts/
cp Output/Aromatic/nonat_DNA_nonat_contacts/2C62.pdb .
cp Output/Aromatic/nonat_DNA_nonat_contacts/polyT.pdb .
cp Output/Aromatic/nonat_DNA_nonat_contacts/run.sh .
mv Output/Aromatic/nonat_DNA_nonat_contacts/polyT.pdb .
mv Output/Aromatic/nonat_DNA_nonat_contacts/2C62_DNA.pdb .
fi

echo "======================================================="
echo "nat_DNA_nonat_contacts: 2C62.pdb as input pdb. Using native DNA sequecne structure as custom input (no notive contacts) along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 2C62_DNA.pdb"
cp ../2C62_DNA.pdb .
echo "Ready to proceed [yes=1/no=0]? (0)"
#read input
input=0
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 2C62_DNA.pdb

mkdir Output/Aromatic/nat_DNA_nonat_contacts
rm -r Output/Aromatic/nat_DNA_nonat_contacts/*
mv * Output/Aromatic/nat_DNA_nonat_contacts/
cp Output/Aromatic/nat_DNA_nonat_contacts/2C62.pdb .
cp Output/Aromatic/nat_DNA_nonat_contacts/2C62_DNA.pdb .
cp Output/Aromatic/nat_DNA_nonat_contacts/run.sh .
mv Output/Aromatic/nat_DNA_nonat_contacts/polyT.pdb .
mv Output/Aromatic/nat_DNA_nonat_contacts/2C62_DNA.pdb .
fi

echo "No Cross terms"
echo "======================================================="
echo "nat_DNA_nat_contacts: 2C62.pdb as input pdb. Using different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface"

echo "Ready to proceed [yes=1/no=0]? (0)"
##read input
input=0i
nput=1
if [ $input -eq 1 ]
then
bash ../../generate_NucproInterface_inputfile.sh
echo "
False
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface #--custom_nuc polyT.pdb

mkdir Output/Blank_control/nat_DNA_nat_contacts
rm -r Output/Blank_control/nat_DNA_nat_contacts/*
mv * Output/Blank_control/nat_DNA_nat_contacts/
cp Output/Blank_control/nat_DNA_nat_contacts/2C62.pdb .
cp Output/Blank_control/nat_DNA_nat_contacts/run.sh .
mv Output/Blank_control/nat_DNA_nat_contacts/polyT.pdb .
mv Output/Blank_control/nat_DNA_nat_contacts/2C62_DNA.pdb .
fi

echo "======================================================="
echo "nonat_DNA_nonat_contacts: 2C62.pdb as input pdb. Using different DNA structure along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "The DNA sequence used is PolyT of lenght 11bases (same as the native bound DNA strand)"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb"
cp ../polyT.pdb .
echo "Ready to proceed [yes=1/no=0]? (0)"
##read input
input=0i
nput=1
if [ $input -eq 1 ]
then
bash ../../generate_NucproInterface_inputfile.sh
echo "
False
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc polyT.pdb

mkdir Output/Blank_control/nonat_DNA_nonat_contacts
rm -r Output/Blank_control/nonat_DNA_nonat_contacts/*
mv * Output/Blank_control/nonat_DNA_nonat_contacts/
cp Output/Blank_control/nonat_DNA_nonat_contacts/2C62.pdb .
cp Output/Blank_control/nonat_DNA_nonat_contacts/polyT.pdb .
cp Output/Blank_control/nonat_DNA_nonat_contacts/run.sh .
mv Output/Blank_control/nonat_DNA_nonat_contacts/polyT.pdb .
mv Output/Blank_control/nonat_DNA_nonat_contacts/2C62_DNA.pdb .
fi

echo "======================================================="
echo "nat_DNA_nonat_contacts: 2C62.pdb as input pdb. Using native DNA sequecne structure as custom input (no notive contacts) along with different parameters for Protein-DNA interface."
echo "In the case an interface file is used: nucpro_interface.input. Thefile is used for providing custom repulsion term when accounting for DNA-Protein interaction"
echo "Command used: $ python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 2C62_DNA.pdb"
cp ../2C62_DNA.pdb .
echo "Ready to proceed [yes=1/no=0]? (0)"
##read input
input=0i
nput=1
if [ $input -eq 1 ]
then
bash ../../generate_NucproInterface_inputfile.sh
echo "
False
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

python ../../nucproSBM.py --aa_pd 2C62.pdb --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --no_Pcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --all_chain --interface --custom_nuc 2C62_DNA.pdb

mkdir Output/Blank_control/nat_DNA_nonat_contacts
rm -r Output/Blank_control/nat_DNA_nonat_contacts/*
mv * Output/Blank_control/nat_DNA_nonat_contacts/
cp Output/Blank_control/nat_DNA_nonat_contacts/2C62.pdb .
cp Output/Blank_control/nat_DNA_nonat_contacts/2C62_DNA.pdb .
cp Output/Blank_control/nat_DNA_nonat_contacts/run.sh .
mv Output/Blank_control/nat_DNA_nonat_contacts/polyT.pdb .
mv Output/Blank_control/nat_DNA_nonat_contacts/2C62_DNA.pdb .
fi
