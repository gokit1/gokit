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

* Notes on eSBM:
Find position of C-beta based on type of SBM used. BT model keeps glycines intact. 
Add Statistical potentials Miyazawa-Jernighan and Betancourt-Thirumalai. CB radii are automatically assigned. Interaction type and strength is automatically assigned. Cut-off based contact-map is generated. Bonds, angles, dihedrals and  distances are added. File conversion to .gro, .top and SBM.INP.



### Installing stuff
Make sure pip is installed and running.
```
$ git clone https://github.com/gokit1/gokit.git
$ cd gokit
$ ./INSTALL
```

## Examples
### Generating a contact-mao

### Flavours





---
GO-kit is mainly developed on a 64-bit OS X machine. 
