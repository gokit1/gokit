echo "The script will generate test files for 2es2."
echo "3 folders will be created 1) control (native complex), 2) native_nuc_test and 3) non_native_nuc_test"
echo "Press Enter 1 to proceed. Any other input will exit the script"
read input
#input=1
if [ $input -eq 1 ]
then
python ../../nucproSBM.py --aa_pdb 2es2.pdb --P_nKd 0.3 --nKd 0.5 --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --interface --control
mkdir Output/
rm -r Output/*
mv * Output/
cp Output/2es2.pdb .
cp Output/input_DNA.pdb .
cp Output/run.sh .
mv Output control
###############################
python ../../nucproSBM.py --aa_pdb 2es2.pdb --P_nKd 0.3 --nKd 0.5 --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --interface
mkdir Output/
rm -r Output/*
mv * Output/
cp Output/2es2.pdb .
cp Output/input_DNA.pdb .
cp Output/run.sh .
cp -r Output/control .
mv Output native_nuc_test

#########
python ../../nucproSBM.py --aa_pdb 2es2.pdb --P_nKd 0.3 --nKd 0.5 --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --interface --custom_nuc input_DNA.pdb
mkdir Output/
rm -r Output/*
mv * Output/
cp Output/2es2.pdb .
cp Output/input_DNA.pdb .
cp Output/run.sh .
cp -r Output/control .
cp -r Output/native_nuc_test .
mv Output non_native_nuc_test
#=============================

fi
exit
