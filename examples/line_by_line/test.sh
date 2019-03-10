#debug from here
# python conmaps.py --get_pdb 1ris > 1.out
 python gokit.py --w_native 1ris.pdb --skip_glycine >2.out
#one -bead
 python gokit.py --w_native 1ris.pdb --skip_glycine > 3.out
 python gokit.py -attype 1 -aa_pdb 1ris.pdb -skip_glycine > 4.out
 python gokit.py --attype 1 -aa_pdb 1ris.pdb -hphobic -skip_glycine > 5.out
 python gokit.py -attype 1 -aa_pdb 1ris.pdb -dsb -skip_glycine > 6.out 
 python gokit.py -attype 1 -aa_pdb 1ris.pdb -dsb -hphobic -skip_glycine > 7.ouy

# two-bead
 python gokit.py --attype 2 --aa_pdb 1ris.pdb --skip_glycine > 8.out
 python gokit.py --attype 2 --aa_pdb 1ris.pdb -btmap -skip_glycine > 9.out
 python gokit.py --attype 2 --aa_pdb 1ris.pdb -mjmap -skip_glycine > 10.out
 python gokit.py --attype 2 --aa_pdb 1ris.pdb -skip_glycine -CA_rad 3.8 -interactions  > 11.out 
python gokit.py --attype 2 --aa_pdb 1ris.pdb --pl_map --CAcom --Ka 200 --Kb 1 --Kd 40 --skip_glycine --interaction  --CA_rad 4.0 --CA_sep 4 --CB_sep 3 --CAB_sep 3 > 12.out
