#debug from here
# python conmaps.py --get_pdb 1ris
#one -bead
# python gokit.py --w_native 1ris.pdb --skip_glycine
# python gokit.py -attype 1 -aa_pdb 1ris.pdb -skip_glycine
# python gokit.py --attype 1 -aa_pdb 1ris.pdb -hphobic -skip_glycine
# python gokit.py -attype 1 -aa_pdb 1ris.pdb -dsb -skip_glycine
# python gokit.py -attype 1 -aa_pdb 1ris.pdb -dsb -hphobic -skip_glycine
# python gokit.py -attype 1 -aa_pdb 1ris.pdb -dsb -hphobic -skip_glycine

# two-bead
python gokit.py --attype 2 --aa_pdb 1ris.pdb --skip_glycine
python gokit.py --attype 2 --aa_pdb 1ris.pdb -btmap -skip_glycine
python gokit.py --attype 2 --aa_pdb 1ris.pdb -mjmap -skip_glycine
python gokit.py --attype 2 --aa_pdb 1ris.pdb -skip_glycine -CA_rad 3.8 -interactions 

