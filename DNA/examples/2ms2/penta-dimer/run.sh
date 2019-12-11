echo "Press Enter 1 to proceed. Any other input will exit the script"
read input
#input=1
if [ $input -eq 1 ]
then
python ../../../nucproSBM.py --aa_pdb penta-dimer.pdb --P_nKd 0.9 --nKd 1.5 --grotop gromacs.top -pdbgro gromacs.gro  --CA_rad 1.9 --skip_glycine  --CBfar --CBcharge --dielec 70 --iconc 0.01 --CB_rad 1.5 --debye --interface
mkdir Output/
rm -r Output/*
mv * Output/
cp Output/penta-dimer.pdb .
cp Output/input_RNA.pdb .
cp Output/run.sh .
###############################
fi