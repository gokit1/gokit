# Go-kit: Enabling energy landscape explorations of proteins

A set of python tools that assist setup and post-hoc analysis of simulations of proteins with various flavours of Structure-Based Models (SBMs) for molecular dynamics and discrete path sampling schemes. 

### What flavours are included:
* Standard cut-off based CA 
* Standarg cut-off based CA (with desolvation barrier approximation for coarse-grained water molecules)
* Standard cut-off based CA + hydrophobic contacts
* Standard cut-off based CA + hydrophobic contacts + desolvation barrier 
* A basic Cheung-Thirumalai two-bead (CA-CB) coarse-grained model(with angles and dihedrals).
* SOP-SC two-bead coarse-grained model with statistical potentials (Betancourt-Thirumalai) and statistical radii. [here](http://biotheory.umd.edu/Supplementary%20Materials/btmap.dat)
* Miyazawa-Jernighan parameters for two-bead models.(Table 5 [here](https://www.ncbi.nlm.nih.gov/pubmed/10336383))
* Custom interaction potentials for two-bead models. (Define your own residue-residue interactions for two-beaad models) and place CB either at COM or at CB or at farthest atom in side-chain.

### What does Go-kit do?
Starting from a PDBID, it generates input files for both GROMACS and OPTIM/PATHSAMPLE potential files with the SBM of choice. After the runs are complete, it can be used to analyse results as well. 


### Installing stuff
Make sure pip is installed and running.
```
$ git clone https://github.com/gokit1/gokit.git
$ cd gokit
$ chmod 755 INSTALL
$ ./INSTALL
```

## Examples: Input files for GROMACS and OPTIM.
### Generating a contact-map
Flags work with both - and --

```
$ python conmaps.py --get_pdb 1ris
```
Remove HETATM records and ensure residue numbering starts from 1. Ensure all fields in file have space separators.
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
$ python gokit.py -attype 1 -aa_pdb 1ris.pdb -skip_glycine
```
Remove --skip_glycine if Hydrogen atoms are present in the PDB file. CB beads are placed on the Glycine Hydrogen atom. This is done in the two-bead model of Thirumalai for e.g. 

One-bead C-alpha: 12-10 LJ + hydrophobic
```
$ python gokit.py --attype 1 -aa_pdb 1ris.pdb -hphobic -skip_glycine
```
One-bead C-alpha: dsb
```
$ python gokit.py -attype 1 -aa_pdb 1ris.pdb -dsb -skip_glycine
```
This is the desolvation barrier potential of Chan et al. 

One-bead C-alpha: dsb + hydrophobic

```
$ python gokit.py -attype 1 -aa_pdb 1ris.pdb -dsb -hphobic -skip_glycine
```
Note: The dsb+hp potential is not implemented in OPTIM. Use only for gromacs4 runs. dsb works in OPTIM though. 

Two-bead model: Cheung-Thirumalai

```
$ python gokit.py --attype 2 --aa_pdb 1ris.pdb --skip_glycine
```
Two-bead model: Betancourt-Thirumalai

```
$ python gokit.py --attype 2 --aa_pdb 1ris.pdb -btmap -skip_glycine
```
Two-bead model: Miyazawa-Jernighan 

```
$ python gokit.py --attype 2 --aa_pdb 1ris.pdb -mjmap -skip_glycine
```

Two-bead model: Customised side-chain interactions (Beta)

```
$ python gokit.py --attype 2 --aa_pdb 1ris.pdb -skip_glycine -CA_rad 3.8 -interactions 
```

Two-bead model: Customised

```
python gokit.py --attype 2 --aa_pdb 1ris.pdb --pl_map --CAcom --Ka 200 --Kb 1 --Kd 40 --skip_glycine --interactions --CA_rad 4.0 --CA_sep 4 --CB_sep 3 --CAB_sep 3 
```


include file called interactions.dat in the format of mjmap.dat or btmap.dat.

### See examples/line_by_line for test output of each command and flag mentioned in the manuscript

Two folders called MD and PATH are generated. MD contains gromacs.top and gromacs.gro file that can be used directly for MD runs. 

See [OPTIM](http://www-wales.ch.cam.ac.uk/OPTIM.doc/node1.html) and [PATHSAMPLE](https://wikis.ch.cam.ac.uk/ro-walesdocs/wiki/index.php/PATHSAMPLE) documentation for generating disconnectivity graphs.

[Example](http://www-wales.ch.cam.ac.uk/examples/OPTIM/t3/) of an OPTIM run for S6 protein.
See [GROMACS](http://www.gromacs.org/Documentation/Installation_Instructions_4.5) documentation. 



Reference: [Go-Kit: A Tool To Enable Energy Landscape Exploration of Proteins](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00007)
--------
GO-kit is mainly developed and tested on a 64-bit MacOS(MOJAVE) machine. GROMACS v4.5.3 and PATHSAMPLE/OPTIM versions on 1.3.2019 were used.
