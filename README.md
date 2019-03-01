# Go-kit: Enabling energy landscape explorations of proteins

A set of python tools that assist setup and post-hoc analysis of simulations of proteins with various flavours of Structure-Based Models (SBMs) for GROMACS and OPTIM runs. 

### What flavours are included:
* Standard cut-off based CA 
* Standard cut-off based CA + hydrophobic contacts
* Standard cut-off based CA + hydrophobic contacts + desolvation barrier 
* A basic Cheung-Thirumalai two-bead (CA-CB) coarse-grained model(with angles and dihedrals).
* SOP-SC two-bead coarse-grained model with statistical potentials (Betancourt-Thirumalai) and statistical radii.:
* Miyazawa-Jernighan parameters for two-bead models.



### What does Go-kit do?
Starting from a PDBID, it generates input files for both GROMACS and OPTIM/PATHSAMPLE potential files with the SBM of choice. After the runs are complete, it can be used to analyse results as well. 


### Installing stuff
Make sure pip is installed and running.
```
$ git clone https://github.com/gokit1/gokit.git
$ cd gokit
$ ./INSTALL
```

## Examples: Input files for GROMACS and OPTIM.
### Generating a contact-map
Flags work with both - and --

```
$ python conmaps.py --get_pdb 1ris.pdb
```
Remove HETATM records and ensure residue numbering starts from 1. 
```
$ python conmaps.py --gconmap 1ris.pdb
```

### Flavours
Generate coarse-grained representation
```
$ python gokit.py --w_native 1ris.pdb --skip_glycine
```
One-bead C-alpha: 12-10 LJ
```
$ python gokit.py -attype 1 -aa_pdb 1ris.pdb 
```
Remove --skip_glycine if Hydrogen atoms are present in the PDB file. CB beads are placed on the Glycine Hydrogen atom. This is done in the two-bead model of Thirumalai for e.g. 

One-bead C-alpha: 12-10 LJ + hydrophobic
```
$ python gokit.py --attype 1 -aa_pdb 1ris.pdb -hphobic 
```
One-bead C-alpha: dsb
```
$ python gokit.py -attype 1 -aa_pdb 1ris.pdb -dsb
```
This is the desolvation barrier potential of Chan et al. 

One-bead C-alpha: dsb + hydrophobic

```
$ python gokit.py -attype 1 -aa_pdb 1ris.pdb -dsb -hphobic
```
Note: The dsb+hp potential is not implemented in OPTIM. Use only for gromacs4 runs. dsb works in OPTIM though. 

Two-bead model: Cheung-Thirumalai. 
```
$ python gokit.py --attype 2 --aa pdb 1ris.pdb --skip_glycine
```
Two-bead model: Betancourt-Thirumalai
```
$ python gokit.py --attype 2 --aa pdb 1ris.pdb -btmap -skip_glycine
```
Two-bead model: Miyazawa-Jernighan 
```
$ python gokit.py --attype 2 --aa pdb 1ris.pdb -mjmap -skip_glycine
```


Two folders called MD and PATH are generated. MD contains gromacs.top and gromacs.gro file that can be used directly for MD runs. 
PATH contains an SBM.INP file that is required for OPTIM runs.

See [OPTIM](http://www-wales.ch.cam.ac.uk/OPTIM.doc/node1.html) and [PATHSAMPLE](https://wikis.ch.cam.ac.uk/ro-walesdocs/wiki/index.php/PATHSAMPLE)

[Example](http://www-wales.ch.cam.ac.uk/examples/OPTIM/t3/) of a sample OPTIM run for S6 that uses an SBM.INP generated form Go-kit. 





--------
GO-kit is mainly developed on a 64-bit OS X machine. 
