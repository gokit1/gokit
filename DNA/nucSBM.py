#!/usr/bin/env python

"""
PyeSBM:
usage: python nucSBM.py --options

"""
from __future__ import print_function
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO,Select
import collections
import numpy as np
from protSBM import protsbm
from util import Utils

class PrePDB(Select):
	def __init__(self):
		return
	def accept_chain(self, chain):
		#select DNA/RNA chain
		if chain in chain_list:
			return 1
		else:
			return 0
	def readPDB(self,file):
		#read pdb as list of lines
		fin = open(file)
		all = fin.read()
		pdblines = all.split("\n")
		pdblines = [x.rstrip().lstrip() for x in pdblines if x.startswith("ATOM")]
		return pdblines
	def terminalResidues(self,file):
		#return terminal residue number
		fin = open(file); all=fin.readlines(); fin.close()
		resnum = 0; terminal_residues = list()
		for i in all:
			if i.startswith("ATOM"):
				resnum = int(i[22:26])
			elif i.startswith("TER"):
				terminal_residues.append(resnum)
		return terminal_residues
	def homoNmer2Monomer(self,file):
		#Convert a homo-multimer to a monomer
		print (">>> converting homo-N-mer(if the pdb is) to a monomer.")
		#reading pdb
		P = PDBParser(PERMISSIVE=0)
		print (file)
		structure = P.get_structure('test', file)
		model = structure.get_chains()
		in_chain = dict()
		residues_in_chain = list()
		global chain_list; chain_list = list()	
		for chain in model:
			rescount = 0 	#initializing chain residue count
			l = list()
			c = chain.get_id()
			for residue in chain:
				#adding residue name to list l
				if [residue.get_resname(),rescount] not in l:
					l.append([residue.get_resname(),rescount])
				rescount += 1
			if l not in residues_in_chain:
				#adding the list of residues in the chain list
				#storing only unique lists
				residues_in_chain.append(l)
				in_chain[c] = l		
				chain_list.append(chain)
		#writing back to file
		io = PDBIO()
		io.set_structure(structure)
		io.save("mono_"+file,PrePDB())
		return "mono_"+file
	def appendChain(self,file):
		#returns chain with renumbered residues
		print (">>> Appending the residues of the chains.\n>>> The chain_id remains same,", file)
		chain_terminal = list()
		fin = open(file)
		#lines = [l for l in fin.readlines() if not l[12:16].strip().startswith('H')]
		lines = [l for l in fin.readlines()]
		nuc_res = ("A","G","T","U","C","DA","DG","DT","DC")
		nuclines = list(); aalines = list()
		multi_occ = dict()	#pdb lines hashed to atom name
		for l in lines:
			if l.startswith(("ATOM","TER")):
				if l.startswith("ATOM"):
					residue = l[17:21].strip()
					if residue in nuc_res:
						nuclines.append(l)
					else:
						if l[16] == " ": #not A,B,C....
							#no multiple occ
							if len(multi_occ)!=0:
								#if there are multiple occupancy
								atnmae = ""
								for atnmae in multi_occ:
									multi_occ[atnmae].sort()	#sorting
									#loadig the one with maximum occ
									prev_line = multi_occ[atnmae][len(multi_occ[atnmae])-1]
									aalines.append(prev_line[1])
								multi_occ = dict()
							aalines.append(l)
						else:
							atname = l[12:16].strip()
							if atname not in multi_occ:
								multi_occ[atname] = list()	
							multi_occ[atname].append((float(l[54:60]),l))
				elif l.startswith("TER"):
					if residue in nuc_res:
						nuclines.append(l.strip())
					else:
						aalines.append(l.strip())					
		if len(nuclines) == 0:
			print ("The PDB file lacks DNA/RNA. Prefer using the protein obly model......The code will continue")
		else:
			if not nuclines[len(nuclines)-1].startswith("TER"):
				#if nucleotide chain doesnot end with a ter
				nuclines.append("TER")
		lines = nuclines + aalines
		if len(lines) == 0:
			print ("Encountered enpty file!!!!!! ",file," returning 0. The code will run normally if this function was not called for the first time")
			return (0,0)
		rescount= 0; atcount = 0
		chain_number = ord(lines[0][21]) #starting with the chain id of the first residue
		prev_res = 0
		outfile = "concat_"+file
		fout = open(outfile,"w+")
		for l in lines:
			if l.startswith("ATOM"):
				if l[22:26]!=prev_res:
					rescount = rescount + 1
				atcount = atcount + 1
				fout.write(l[0:6]+str(atcount).rjust(5)+l[11:21]+chr(chain_number)+str(rescount).rjust(4)+l[26:])
				prev_res = l[22:26]
			elif l.startswith("TER"):
				fout.write("TER".ljust(6)+str(atcount).rjust(5)+l[11:21].ljust(10)+chr(chain_number)+str(rescount).rjust(4)+"\n")
				chain_terminal.append(rescount)
				chain_number += 1
			else:
				fout.write(l)    
		fout.close()
		return [outfile,chain_terminal]
	def sepNucPro(self,file):
		#Function separates nucleotide and protein components of a pdbfile
		P = PDBParser(PERMISSIVE=0)
		structure = P.get_structure('test', file)
		print (">>> writing nucleotide only file: ","nuc_"+file)
		io = PDBIO()
		io.set_structure(structure)
		io.save("nuc_"+file,nucsbm())
		print (">>> writing protein only file: ","aa_"+file)
		io.set_structure(structure)
		io.save("aa_"+file,protsbm())
		return ("nuc_"+file,"aa_"+file)
	def atomIn(self):
		#dictionary of atoms in coarse grain beads
		atoms_in = dict()
		atoms_in["CA"] = ("N","C","CA","O")
		#atoms_in["CB"]  = ['CB', 'CD', 'CD1', 'CD2', 'CD3', 'CE', 'CE1', 'CE2', 'CE3', 'CG', 'CG1', 'CG2', 'CG3', 'CH', 'CH1','CH2', 'CH', 'CZ', 'CZ1OG', 'CZ2', 'CZ3', 'ND1', 'ND2', 'NE', 'NE1NE2', 'NH1', 'NH2', 'NZ', 'OD1', 'OD2', 'OE1', 'OE2', 'OG1', 'OH', 'SD', 'SG']
		atoms_in["B"] = ("N1","C2","N2","N3","C4","C5","C6","O6","N6","N7","C8","N9","O2","O4","N4","C7")
		atoms_in["S"] = ("C1'","C2'","C3'","C4'","C5'","O2'","O3'","O4'","O5'")
		atoms_in["XP"] = ("P","OP1","OP2","O1P","O2P")
		return atoms_in
	def stackEps(self):
		#Value of epsilon for stacking potential (represented as a LJ term)
		#the values down here are relative to the native contact strength i.e. 1KJ/mol

		#reference for base-base stacking		#https://doi.org/10.1371/journal.pcbi.1006768.t002
		#reference for base-CB stacking			#https://doi.org/10.1371/journal.pcbi.1006768.g001
		d1 =	{
			'BT': {'BG': '1.7', 'BA': '1.4', 'BC': '1.4', 'BT': '1.0', 'BU': '1.0', 'W' : '4.3', 'Y': '2.3', 'F' : '2.8', 'H' : '2.5'},
			'BU': {'BG': '1.5', 'BA': '1.3', 'BC': '1.3', 'BT': '1.0', 'BU': '0.9', 'W' : '3.9', 'Y': '2.1', 'F' : '2.6', 'H' : '2.3'},
			'BG': {'BG': '1.8', 'BA': '1.8', 'BC': '1.5', 'BT': '1.7', 'BU': '1.5', 'W' : '3.9', 'Y': '2.9', 'F' : '3.0', 'H' : '2.2'},
			'BA': {'BG': '1.8', 'BA': '1.4', 'BC': '1.5', 'BT': '1.4', 'BU': '1.3', 'W' : '3.1', 'Y': '3.3', 'F' : '3.3', 'H' : '2.2'}, 
			'BC': {'BG': '1.5', 'BA': '1.5', 'BC': '1.4', 'BT': '1.4', 'BU': '1.3', 'W' : '2.4', 'Y': '1.4', 'F' : '1.7', 'H' : '1.4'}}
		return d1
		#d = dict()
		#f=open("stackEpsilon.dat")
		#for i in f:
		#	if i[0] in ("A","T","G","C","U"):
		#		a1 = "B" + i[0]
		#		if a1 not in d:
		#			d[a1] = dict()
		#		a2=i[1]
		#		if a2 in ("A","T","G","C","U"):
		#			a2 = "B" + a2
		#		d[a1][a2] = float(i.split()[1])
		#return d
	def basepairEps(self):
		#defining LJ strenfth for base pairs (will be effectively 4 times of the below value)
		#theese values are  parameterized over multiple runs
		return {("A","T"):1.0, ("T","A"):1.0, ("A","U"):1.0, ("U","A"):1.0, ("G","C"):1.5,("C","G"):1.5}

class Measure():
	#__________________TO___BE____SHIFTED________#
	def distances(self,file,pairs):
		#loading co-ordinates
		fin = open(file)
		xyz = list()
		for i in fin:
			if i.startswith("ATOM"):
				xyz.append(np.float_([i[30:38],i[38:46],i[46:54]]))
		xyz = np.array(xyz)
		fin.close()

		#calcularing bond lengths
		r = list()
		for i in pairs:
			r.append(sum((xyz[i[1]]-xyz[i[0]])**2)**0.5)
		
		#storing bond lengths along with bond indexes
		dist = list()
		dist.append(r)
		for i in pairs:
			dist.append(i)
		
		return dist
	def angles(self,file,triplets):
		#loading co-ordinates
		fin = open(file)
		xyz = list()
		for i in fin:
			if i.startswith("ATOM"):
				xyz.append(np.float_([i[30:38],i[38:46],i[46:54]]))
		xyz = np.array(xyz)
		fin.close()
		
		#calculating angles
		theta = list()
		for t in triplets:
			#angle xzy0,xyz1,xyz2
			AB = xyz[t[0]]-xyz[t[1]]; BC = xyz[t[2]]-xyz[t[1]]
			AB_BC = (sum(AB**2)**0.5)*(sum(BC**2)**0.5)
			theta.append(np.arccos(sum(AB*BC)/AB_BC))	#	in rads
		
		#returning angles along with teir triplets
		angles = list()
		angles.append(theta)
		for i in triplets:
			angles.append(i)
		return angles
	def dihedrals(self,file,quads):
		#loading co-ordinates
		fin = open(file)
		xyz = list()
		for i in fin:
			if i.startswith("ATOM"):
				xyz.append(np.float_([i[30:38],i[38:46],i[46:54]]))
		xyz = np.array(xyz)
		fin.close()

		#calculating torsions
		phi = list()
		for q in quads:
			#vectors
			BA = xyz[q[0]]-xyz[q[1]]; BC = xyz[q[2]]-xyz[q[1]]
			CD = xyz[q[3]]-xyz[q[2]]; CB = BC

			#normals to plane 1 and 2
			#n1 = np.array([BA[1]*BC[2]-BA[2]*BC[1],BA[2]*BC[0]-BA[0]*BC[2],BA[0]*BC[1]-BA[1]*BC[0]])
			#n2 = np.array([CD[1]*CB[2]-CD[2]*CB[1],CD[2]*CB[0]-CD[0]*CB[2],CD[0]*CB[1]-CD[1]*CB[0]])
			n1 = np.cross(BA,BC)
			n2 = np.cross(CD,BC)


			#checking the sign (+/-) of the angle
				#Cos-1(x) usually has 2 solutions (+phi) and (-phi)
				#By default acos (from math lib) or arccos (from numpy) 
				#will return only +phi as they are calculated in the limits
				# if 0<phi<pi. This will give a mirrored dihedrals
				#for places where phi was actually negative
			#a vector prependicular to both n1 and m2 will be
			#parallel or anti-parallel to BC

			direction = np.cross(n1,n2)
			cosine_d = sum(direction*BC)/((sum(direction**2)**0.5)*(sum(BC**2)**0.5))

			#if cosine_d is -1, the n1Xn2 is antiparallel to BC and the dihedral is negative
			#if cosine_d is  1, the n1Xn2 is parallel to BC and the dihedral is positive		 
			
			sign = int(np.round(cosine_d))
			#making dihedtals are negative, they lie b/w pi to 2*pi positive
			if sign == 1:
				angle_normalizer = 0
			elif sign == -1:
				angle_normalizer = 2*np.pi

			#calculating dihedral
			n1_n2 = (sum(n1**2)**0.5)*(sum(n2**2)**0.5)
			phi.append(angle_normalizer + sign*np.arccos(sum(n1*n2)/n1_n2))	#	in rads

		#returning dihedrals along with teir quadruplets
		dihedrals = list()
		dihedrals.append(phi)
		for i in quads:
			dihedrals.append(i)
		return dihedrals

class nucsbm(protsbm,Select):
	def __init__(self):
		return
	def accept_residue(self,residue):
		#accept only nucleotide residues
		global nuc_res 
		nuc_res = ("A","G","T","U","C","DA","DG","DT","DC")
		if residue.get_resname().lstrip().rstrip() in nuc_res:
			return 1
		else:
			return 0
	def writeSample_BDNA(self,pdbfile):
		#sample B-form dsDNA (Genaretd from NUCGEN+. Please refer to Manual for more details)
		pdb = ["ATOM      1  C1'  DA A   1       5.778   1.547   1.713  1.00     0\n", 'ATOM      2  N9   DA A   1       4.567   2.399   1.562  1.00     0\n', 'ATOM      3  C8   DA A   1       4.525   3.728   1.329  1.00     0\n', 'ATOM      4  N7   DA A   1       3.282   4.193   1.246  1.00     0\n', 'ATOM      5  C5   DA A   1       2.503   3.084   1.438  1.00     0\n', 'ATOM      6  C4   DA A   1       3.266   1.978   1.633  1.00     0\n', 'ATOM      7  C6   DA A   1       1.101   2.919   1.464  1.00     0\n', 'ATOM      8  N6   DA A   1       0.214   3.919   1.287  1.00     0\n', 'ATOM      9  N1   DA A   1       0.622   1.664   1.682  1.00     0\n', 'ATOM     10  C2   DA A   1       1.507   0.658   1.860  1.00     0\n', 'ATOM     11  N3   DA A   1       2.853   0.700   1.855  1.00     0\n', "ATOM     12  O4'  DA A   1       6.391   1.654   2.990  1.00     0\n", "ATOM     13  O5'  DA A   1       8.914   3.477   2.623  1.00     0\n", "ATOM     14  C5'  DA A   1       8.679   2.365   3.506  1.00     0\n", "ATOM     15  C4'  DA A   1       7.803   1.335   2.821  1.00     0\n", "ATOM     16  C3'  DA A   1       8.060   1.130   1.328  1.00     0\n", "ATOM     17  O3'  DA A   1       7.985  -0.260   1.030  1.00     0\n", "ATOM     18  C2'  DA A   1       6.868   1.843   0.690  1.00     0\n", 'ATOM     19  P    DA A   1       9.131  -0.927   0.134  1.00     0\n', 'ATOM     20  OA   DA A   1       9.485  -0.015  -0.978  1.00     0\n', 'ATOM     21  OB   DA A   1      10.248  -1.361   1.002  1.00     0\n', "ATOM     42  C1'  DG A   2       4.703  -1.913  -1.572  1.00     0\n", 'ATOM     43  N9   DG A   2       4.227  -0.510  -1.729  1.00     0\n', 'ATOM     44  C8   DG A   2       4.980   0.590  -1.943  1.00     0\n', 'ATOM     45  N7   DG A   2       4.261   1.699  -2.039  1.00     0\n', 'ATOM     46  C5   DG A   2       2.955   1.233  -1.868  1.00     0\n', 'ATOM     47  C6   DG A   2       1.720   1.956  -1.869  1.00     0\n', 'ATOM     48  C4   DG A   2       2.902  -0.093  -1.679  1.00     0\n', 'ATOM     49  O6   DG A   2       1.537   3.156  -2.021  1.00     0\n', 'ATOM     50  N1   DG A   2       0.623   1.122  -1.664  1.00     0\n', 'ATOM     51  C2   DG A   2       0.680  -0.246  -1.478  1.00     0\n', 'ATOM     52  N2   DG A   2      -0.459  -0.909  -1.294  1.00     0\n', 'ATOM     53  N3   DG A   2       1.861  -0.925  -1.479  1.00     0\n', "ATOM     54  O4'  DG A   2       5.236  -2.190  -0.285  1.00     0\n", "ATOM     55  O5'  DG A   2       8.355  -2.199  -0.590  1.00     0\n", "ATOM     56  C5'  DG A   2       7.494  -2.964   0.274  1.00     0\n", "ATOM     57  C4'  DG A   2       6.193  -3.279  -0.438  1.00     0\n", "ATOM     58  C3'  DG A   2       6.310  -3.591  -1.930  1.00     0\n", "ATOM     59  O3'  DG A   2       5.437  -4.670  -2.249  1.00     0\n", "ATOM     60  C2'  DG A   2       5.778  -2.311  -2.575  1.00     0\n", 'ATOM     61  P    DG A   2       5.990  -5.881  -3.138  1.00     0\n', 'ATOM     62  OA   DG A   2       6.834  -5.348  -4.231  1.00     0\n', 'ATOM     63  OB   DG A   2       6.620  -6.892  -2.261  1.00     0\n', 'TER\n', 'ATOM     82  OB   DC B   3     -10.069  -1.635  -0.950  1.00     0\n', 'ATOM     81  OA   DC B   3      -9.325  -0.162   0.945  1.00     0\n', 'ATOM     80  P    DC B   3      -8.958  -1.133  -0.111  1.00     0\n', "ATOM     79  C2'  DC B   3      -6.741   1.635  -0.835  1.00     0\n", "ATOM     78  O3'  DC B   3      -7.824  -0.501  -1.048  1.00     0\n", "ATOM     77  C3'  DC B   3      -7.922   0.868  -1.427  1.00     0\n", "ATOM     76  C4'  DC B   3      -7.670   0.988  -2.931  1.00     0\n", "ATOM     75  C5'  DC B   3      -8.563   1.962  -3.673  1.00     0\n", "ATOM     74  O5'  DC B   3      -8.814   3.120  -2.856  1.00     0\n", "ATOM     73  O4'  DC B   3      -6.264   1.319  -3.121  1.00     0\n", 'ATOM     72  O2   DC B   3      -3.009   0.340  -1.956  1.00     0\n', 'ATOM     71  C2   DC B   3      -3.066   1.567  -1.815  1.00     0\n', 'ATOM     70  N3   DC B   3      -1.969   2.340  -1.728  1.00     0\n', 'ATOM     69  N4   DC B   3      -0.953   4.391  -1.495  1.00     0\n', 'ATOM     68  C4   DC B   3      -2.049   3.647  -1.579  1.00     0\n', 'ATOM     67  C5   DC B   3      -3.337   4.298  -1.502  1.00     0\n', 'ATOM     66  C6   DC B   3      -4.429   3.517  -1.589  1.00     0\n', 'ATOM     65  N1   DC B   3      -4.449   2.176  -1.743  1.00     0\n', "ATOM     64  C1'  DC B   3      -5.647   1.297  -1.841  1.00     0\n", 'ATOM     41  OB   DT B   4      -6.319  -6.909   2.631  1.00     0\n', 'ATOM     40  OA   DT B   4      -6.556  -5.254   4.506  1.00     0\n', 'ATOM     39  P    DT B   4      -5.704  -5.837   3.445  1.00     0\n', "ATOM     38  C2'  DT B   4      -5.550  -2.305   2.669  1.00     0\n", "ATOM     37  O3'  DT B   4      -5.172  -4.673   2.484  1.00     0\n", "ATOM     36  C3'  DT B   4      -6.062  -3.629   2.103  1.00     0\n", "ATOM     35  C4'  DT B   4      -5.951  -3.405   0.594  1.00     0\n", "ATOM     34  C5'  DT B   4      -7.257  -3.154  -0.133  1.00     0\n", "ATOM     33  O5'  DT B   4      -8.129  -2.353   0.685  1.00     0\n", "ATOM     32  O4'  DT B   4      -5.011  -2.313   0.375  1.00     0\n", 'ATOM     31  O2   DT B   4      -1.762  -1.142   1.478  1.00     0\n', 'ATOM     30  C2   DT B   4      -2.588  -0.254   1.615  1.00     0\n', 'ATOM     29  N3   DT B   4      -2.172   1.060   1.682  1.00     0\n', 'ATOM     28  O4   DT B   4      -2.524   3.296   1.884  1.00     0\n', 'ATOM     27  C4   DT B   4      -3.027   2.122   1.835  1.00     0\n', 'ATOM     26  C7   DT B   4      -5.449   2.931   2.099  1.00     0\n', 'ATOM     25  C5   DT B   4      -4.417   1.843   1.929  1.00     0\n', 'ATOM     24  C6   DT B   4      -4.799   0.565   1.862  1.00     0\n', 'ATOM     23  N1   DT B   4      -4.029  -0.534   1.714  1.00     0\n', "ATOM     22  C1'  DT B   4      -4.482  -1.950   1.642  1.00     0\n", 'END']
		fout = open(pdbfile,"w+")
		for i in pdb:
			fout.write(i)
		fout.close()
		return 1
	def writeSample_ARNA(self,pdbfile):
		pdb = ["ATOM      1  C1'   A A   1      -6.296  -0.042   2.314  1.00     0\n", 'ATOM      2  N9    A A   1      -5.206  -0.987   2.100  1.00     0\n', 'ATOM      3  C8    A A   1      -5.296  -2.349   1.937  1.00     0\n', 'ATOM      4  N7    A A   1      -4.138  -2.937   1.761  1.00     0\n', 'ATOM      5  C5    A A   1      -3.223  -1.894   1.812  1.00     0\n', 'ATOM      6  C4    A A   1      -3.867  -0.688   2.020  1.00     0\n', 'ATOM      7  C6    A A   1      -1.823  -1.859   1.694  1.00     0\n', 'ATOM      8  N6    A A   1      -1.068  -2.941   1.492  1.00     0\n', 'ATOM      9  N1    A A   1      -1.216  -0.656   1.792  1.00     0\n', 'ATOM     10  C2    A A   1      -1.973   0.428   1.995  1.00     0\n', 'ATOM     11  N3    A A   1      -3.295   0.524   2.123  1.00     0\n', "ATOM     12  O4'   A A   1      -7.410  -0.735   2.856  1.00     0\n", "ATOM     13  O5'   A A   1      -8.489  -3.072   1.353  1.00     0\n", "ATOM     14  C5'   A A   1      -9.255  -2.025   1.982  1.00     0\n", "ATOM     15  C4'   A A   1      -8.506  -0.711   1.895  1.00     0\n", "ATOM     16  C3'   A A   1      -7.820  -0.423   0.559  1.00     0\n", "ATOM     17  O3'   A A   1      -8.714   0.161  -0.378  1.00     0\n", "ATOM     18  C2'   A A   1      -6.800   0.614   1.032  1.00     0\n", "ATOM     19  O2'   A A   1      -7.226   1.934   1.305  1.00     0\n", 'ATOM     20  P     A A   1      -9.024  -4.578   1.432  1.00     0\n', 'ATOM     21  OA    A A   1     -10.312  -4.686   0.713  1.00     0\n', 'ATOM     22  OB    A A   1      -7.956  -5.506   0.997  1.00     0\n', "ATOM     43  C1'   G A   2      -3.826   2.982  -1.443  1.00     0\n", 'ATOM     44  N9    G A   2      -3.381   1.592  -1.442  1.00     0\n', 'ATOM     45  C8    G A   2      -4.169   0.470  -1.346  1.00     0\n', 'ATOM     46  N7    G A   2      -3.483  -0.640  -1.373  1.00     0\n', 'ATOM     47  C5    G A   2      -2.163  -0.228  -1.495  1.00     0\n', 'ATOM     48  C6    G A   2      -0.967  -0.985  -1.575  1.00     0\n', 'ATOM     49  C4    G A   2      -2.084   1.148  -1.539  1.00     0\n', 'ATOM     50  O6    G A   2      -0.830  -2.215  -1.553  1.00     0\n', 'ATOM     51  N1    G A   2       0.152  -0.166  -1.692  1.00     0\n', 'ATOM     52  C2    G A   2       0.122   1.207  -1.727  1.00     0\n', 'ATOM     53  N2    G A   2       1.310   1.820  -1.842  1.00     0\n', 'ATOM     54  N3    G A   2      -0.986   1.924  -1.654  1.00     0\n', "ATOM     55  O4'   G A   2      -5.094   3.062  -0.811  1.00     0\n", "ATOM     56  O5'   G A   2      -7.312   1.379  -1.873  1.00     0\n", "ATOM     57  C5'   G A   2      -7.397   2.761  -1.477  1.00     0\n", "ATOM     58  C4'   G A   2      -6.094   3.470  -1.790  1.00     0\n", "ATOM     59  C3'   G A   2      -5.456   3.126  -3.136  1.00     0\n", "ATOM     60  O3'   G A   2      -5.949   3.892  -4.195  1.00     0\n", "ATOM     61  C2'   G A   2      -4.022   3.563  -2.840  1.00     0\n", "ATOM     62  O2'   G A   2      -3.700   4.940  -2.835  1.00     0\n", 'ATOM     63  P     G A   2      -8.526   0.396  -1.525  1.00     0\n', 'ATOM     64  OA    G A   2      -9.741   0.833  -2.248  1.00     0\n', 'ATOM     65  OB    G A   2      -8.109  -1.007  -1.746  1.00     0\n', 'TER\n', 'ATOM     85  OB    C B   3       7.827  -5.349  -1.336  1.00     0\n', 'ATOM     84  OA    C B   3      10.194  -4.573  -1.017  1.00     0\n', 'ATOM     83  P     C B   3       8.903  -4.407  -1.720  1.00     0\n', "ATOM     82  O2'   C B   3       7.179   2.104  -1.178  1.00     0\n", "ATOM     81  C2'   C B   3       6.740   0.773  -0.986  1.00     0\n", "ATOM     80  O3'   C B   3       8.607   0.261   0.380  1.00     0\n", "ATOM     79  C3'   C B   3       7.750  -0.300  -0.584  1.00     0\n", "ATOM     78  C4'   C B   3       8.426  -0.513  -1.939  1.00     0\n", "ATOM     77  C5'   C B   3       9.160  -1.829  -2.112  1.00     0\n", "ATOM     76  O5'   C B   3       8.386  -2.903  -1.545  1.00     0\n", "ATOM     75  O4'   C B   3       7.325  -0.465  -2.893  1.00     0\n", 'ATOM     74  O2    C B   3       3.649   0.989  -2.109  1.00     0\n', 'ATOM     73  C2    C B   3       3.818  -0.239  -2.050  1.00     0\n', 'ATOM     72  N3    C B   3       2.780  -1.094  -1.899  1.00     0\n', 'ATOM     71  N4    C B   3       1.951  -3.216  -1.689  1.00     0\n', 'ATOM     70  C4    C B   3       3.004  -2.409  -1.839  1.00     0\n', 'ATOM     69  C5    C B   3       4.317  -2.954  -1.928  1.00     0\n', 'ATOM     68  C6    C B   3       5.334  -2.096  -2.077  1.00     0\n', 'ATOM     67  N1    C B   3       5.115  -0.748  -2.140  1.00     0\n', "ATOM     66  C1'   C B   3       6.223   0.205  -2.303  1.00     0\n", 'ATOM     42  OB    U B   4       7.970  -1.174   1.792  1.00     0\n', 'ATOM     41  OA    U B   4       9.609   0.621   2.418  1.00     0\n', 'ATOM     40  P     U B   4       8.396   0.237   1.662  1.00     0\n', "ATOM     39  O2'   U B   4       3.588   4.720   3.223  1.00     0\n", "ATOM     38  C2'   U B   4       3.902   3.344   3.142  1.00     0\n", "ATOM     37  O3'   U B   4       5.856   3.605   4.525  1.00     0\n", "ATOM     36  C3'   U B   4       5.331   2.880   3.421  1.00     0\n", "ATOM     35  O4'   U B   4       4.984   2.964   1.094  1.00     0\n", "ATOM     34  C4'   U B   4       5.980   3.303   2.103  1.00     0\n", "ATOM     33  C5'   U B   4       7.281   2.607   1.754  1.00     0\n", "ATOM     32  O5'   U B   4       7.186   1.203   2.063  1.00     0\n", 'ATOM     31  O2    U B   4       1.075   2.094   1.846  1.00     0\n', 'ATOM     30  C2    U B   4       1.904   1.211   1.702  1.00     0\n', 'ATOM     29  N3    U B   4       1.559  -0.114   1.609  1.00     0\n', 'ATOM     28  O4    U B   4       1.953  -2.328   1.378  1.00     0\n', 'ATOM     27  C4    U B   4       2.410  -1.190   1.446  1.00     0\n', 'ATOM     26  C5    U B   4       3.814  -0.858   1.367  1.00     0\n', 'ATOM     25  C6    U B   4       4.162   0.432   1.457  1.00     0\n', 'ATOM     24  N1    U B   4       3.256   1.458   1.620  1.00     0\n', "ATOM     23  C1'   U B   4       3.712   2.853   1.711  1.00     0\n", 'END']
		fout = open(pdbfile,"w+")
		for i in pdb:
			fout.write(i)
		fout.close()
		return 1
	def getIndices(self,pdbfile):
		#get atomindices
		#load pdbfile
		structure = PrePDB().readPDB(pdbfile)
		CG_indices = dict()
		for l in structure:
			atom_num = int(l[6:11].rstrip().lstrip())
			atom_name = l[12:16].rstrip().lstrip()
			res_id = int(l[22:26].rstrip().lstrip())
			X = float(l[30:38].rstrip().lstrip());Y = float(l[38:46].rstrip().lstrip());Z =	float(l[46:54].rstrip().lstrip())	
			xyz = np.array([X,Y,Z])
			#complete structure hased to atomnumber
			CG_indices[(atom_num)] =[res_id,atom_name,xyz]
		terminal_residues = PrePDB().terminalResidues(pdbfile)
		return CG_indices,terminal_residues
	def getCoordinates(self,pdbfile,CG_pos):
		#get co-rdinates for the CG beads by atom in the bead or COM
		print (">>> Generating Nuleic Acids Coarse-Grain Co-ordinates.")
		#loading pdb file
		P = PDBParser(PERMISSIVE=0)
		structure = P.get_structure('test', pdbfile)

		#atome name dictioanry
		#the atomname are defined in a separate function as well but the discroption is different
		atomname=dict()
		atomname["Bpu"] = ("N1","C2","H2-N2","N3","C4","C5","C6","O6-N6","N7","C8","N9")
		atomname["Bpy"] = ("N1","C2","O2","N3","C4","O4-N4","C5","C6")
		atomname["S"] = ("C1'","C2'","C3'","C4'","C5'","H2'-O2'","O3'","O4'","O5'")
		atomname["P"] = ("P","OP1","OP2","O1P","O2P","OA","OB")

		#Nucleotide residues in the PDB
		nuc_res = ("A","G","T","U","C","DA","DG","DT","DC")
		bases = dict()
		bases["Bpy"] = ("T","C","U","DT","DC") #pyrimidine residues
		bases["Bpu"] = ("A","G","DA","DG")		#purines residue
		
		count = 0
		CG_coord = dict()
		for group,pos in CG_pos.items():
			#looping over group name and user defined position for the bead
			#group = Phosphate P, Sugar S, Purines Pu and Pyrimidines Py
			CG_coord[group] = []
			count=0
			for chain in structure.get_chains():
				#looping over chains in the structure
				res_count = 0		#initializing residues in the chain
				c=chain.get_id()	#chain id
				for residue in chain:
					#looping over residues in the chain
					r=residue.get_resname().strip()
					coord = []	#co-ordinates array
					mass = []	#atomic mass array
					if r in nuc_res:
						#looping over predefined list of nucleotide residues
						if group[0]=="B" and r not in bases[group]:
							res_count = res_count + 1
							#Loop over the next base
							continue	#skipping this loop as group not in residue
						if pos=="COM":
							#if user defined position for the bead is group's center of masss
							#print(pos)
							for atom in residue:
								#looping over atoms in the residue
								a=atom.get_name().strip()
								#print (a,residue.id[1])
								if a in atomname[group]:
									#getting co-ordinates and mass
									coord.append(atom.get_coord())
									mass.append(atom._assign_atom_mass())
							#assert len(mass)!=0, "Zero length mass array!"
							coord = np.array(coord)
							mass = np.array(mass)
							assert len(coord)==len(mass)
							mass_sum=sum(mass)
							if len(mass)==0:
								res_count = res_count + 1
								continue
								#the group is absent in this resiude
							#determining Center of mass
							#0: COM co-ordinates, 1:residue count (counted by counter), 2: chain id, residue count  (from pdb), residue name 
							COM = [[np.matmul(([float(atom_mass / mass_sum) for atom_mass in mass]),coord)],[count+res_count],[(c,residue.id[1],r)]]
							#print(COM)
							#COM is a array of 2 elements in each row: 1)list of co-ordinates. 2) residue number
							CG_coord[group]+=COM
							#print ("COM",CG_coord[group])
						else:
							#if position is not Ceneter of mass
							#if pos != COM, pos will be atom-name from all-atompdb
							#try:
							if pos in residue:
							#try:
								#check if the atom is present in residue
								atom = residue[pos]

								#take the co-ordinates as it is		
								coord.append(atom.get_coord())
								mass.append(atom._assign_atom_mass())
								coord = np.array(coord)
								mass = np.array(mass)
								mass_sum=sum(mass)
								COM = [[np.matmul(([float(atom_mass / mass_sum) for atom_mass in mass]),coord)],[count+res_count],[(c,residue.id[1],r)]]
								CG_coord[group]+=COM
							#except:
							#	print (pos,"is missing from the rediue ",count+res_count," Skipping the bead!!!!")
						res_count = res_count + 1
				#counts the residue in chain and adds it to the count of residue in next chain
				count = count + res_count
				#res_count is redundant. Newer code have a refine function which refines residue index before
				#starting the main program. The index  is used later in the code and hence will not be removed till some further update
			#if group == "P":
			#	print (group,CG_coord[group])
			CG_coord[group]=np.array(CG_coord[group]).reshape(int(len(CG_coord[group])/3),3)
		#if not pdbfile.startswith("sample"):
		#exit()
		return CG_coord
	def writeNative(self,pdbfile,CG_pos,grofile):	
		if len(grofile)==0:
			grofile = pdbfile+".gro"
		print (">>> Writing Nucleic Acids native and .gro files.", pdbfile, grofile)

		#output native nucleotide CG file 
		nuc_native_out = 'nuc_native_P-S-B.pdb'
		#checkiong if the current call is for generating sample files
		if grofile == "b.dna.gro":
			#generating sample all atom file for B-form DNA
			self.writeSample_BDNA(pdbfile)
			nuc_native_out = "b.dna.pdb"
		elif grofile == "a.rna.gro":
			#generating sample all atom file for A-form RNA
			self.writeSample_ARNA(pdbfile)
			nuc_native_out = "a.rna.pdb"


		checksize = open(pdbfile)
		lines = checksize.readlines()
		if len(lines) <= 1:
			print ("As mentioned earlier, the ibput file is PROTEIN ONLY.....Exiting now")
			return 0
		else:
			#if the file is not emplty
			terminal_residues = list()
			for i in lines:
				if i.startswith("ATOM"):
					resnum = int(i[22:26])
				elif i.startswith("TER"):
					terminal_residues.append(resnum)

		#Defining Nitrogen Base types
		bases = dict()
		bases["Bpy"] = ("T","C","U","DT","DC")
		bases["Bpu"] = ("A","G","DA","DG")

		#Getting co-ordinates for beads
		CG_coord =dict()
		
		#function returns co-ordintates, residue, chain_id
		CG_coord = self.getCoordinates(pdbfile,CG_pos)
		CG_others = dict()
		
		natoms_CG = 0		#counter for number of atoms
		res_num  = list()	#list of index for the residue
		for group,CG_parameters in CG_coord.items():
			coord = dict(); others = dict()			
			for (x,y,z) in CG_parameters:
				#separating co-rdinates from other parameters and storing at same index
				coord[int(y)] = x.tolist()      #x = [X,Y,Z] co-ordinates
				others[int(y)] = z				#z = [chain_id,res_number,res_name]			
				if int(y) not in res_num:
					res_num.append(int(y))
				#print (y); #print (coord[y])
				natoms_CG = natoms_CG + 1
			#segrigating based on group
			CG_coord[group] = coord
			CG_others[group] = others
		res_num.sort()
		#Writing PDB
		l = {}	#each line for native .pdb file
		k = {}	#each line for coarse grain .gro file
		iter = 0	#iteration variable

		#for .gro file
		fgro = open(grofile,"w+")
		fgro.write('%s\n'%('Grofile generated from Lab16 repo. See https://bitbucket.org/nsridhar/lab16'))#.write('%s %s%s%s\n' % ("CG-MD for",pdbfile.split(".")[0],"\n",str(natoms_CG).rjust(5)))
		fgro.write('%s\n' % str(natoms_CG).rjust(5))
		
		#for native coarse grain pdb file
		fout = open(nuc_native_out, "w+")
		#for i in range(1,len(CG_coord["S"])+1):
		
		for i in res_num:
			#looping over residue index
			for group in ["P","S","Bpu","Bpy"]:
				#looping over group names
				coord = CG_coord[group]							#loading dictionary of co-ordinates
				if i in CG_coord[group]:
				#try:	#The residue may be existing but
					r = CG_others[group][i][2]
																#residue number
					if group[0]=="B" and r not in bases[group]:
						continue	#skipping this loop as group not in residue

					iter = iter + 1								#counting atom number
					l[0] = "ATOM".ljust(6)						#ATOM 1-6,	6char
					l[1] = str(iter).rjust(5)					#atom number 7-11, 5char
					l[2] = "".ljust(1)							#blank at 12, 1char		
					if group[0] == "S":
						if r[0] == "D":							#if bleongs to DNA
							l[3] = str("D"+group[0]).center(4)	#atom name (P,S or B) 13-16, 4char
						else:									#else belongs to RNA
							l[3] = str("R"+group[0]).center(4)	
					elif group[0] == "P":						#phosphate bead name XP
						l[3] = str("X"+group[0]).center(4)
					elif group[0] == "B":						#base bead
						l[3] = str(group[0]+r.strip()[len(r)-1]).center(4)
						#for DNA residue residue name is DA, for RNA it is A. Last resiude of the string
																#bead name
																
					l[4] = "".ljust(1)							#atom indicator/blank	17, 1char
					l[5] = r.ljust(3)							#residue name 18-20, 3char
					l[6] = "".ljust(1)							#blank 21, 1char
					l[7] = CG_others[group][i][0].ljust(1)		#chain 22, 1char
					l[8] = str(CG_others[group][i][1]).rjust(4)	#residue number 23-26, 4char
					l[9] = "".ljust(4)							# insertinal code 27, blank 28-30, 4char
					#print (CG_coord[group][i][0][0])
					l[10]= str('%8.3f' % float(CG_coord[group][i][0])).rjust(8)	#X 31-38, 8char 
					l[11]= str('%8.3f' % float(CG_coord[group][i][1])).rjust(8)	#X 39-46, 8char
					l[12]= str('%8.3f' % float(CG_coord[group][i][2])).rjust(8)	#X 47-54, 8char
					l[13]= "1.00".rjust(6)							#occuplancy factor 55-61, 6char
					l[14]= "00.00".rjust(6)							#temperature factor 62-68, 6char
					l[15]= "".rjust(12)							# element/labels/foontnotes name 69-80
					#print (CG_coord[group][i][0][0])
					#writing grofile
					#print("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n" % (l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14],l[15]))
					fout.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n" % (l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14],l[15]))
					k[0] = l[8].rjust(5)						#residue number grofile
					k[1] = l[5].ljust(5)						#residue name
					k[2] = l[3].rstrip().lstrip().rjust(5)		#atom name
					k[3] = l[1]									#atom number
					k[4]= str('%8.3f' % (float(CG_coord[group][i][0])/10)).rjust(8)	#X 31-38, 8char 
					k[5]= str('%8.3f' % (float(CG_coord[group][i][1])/10)).rjust(8)	#X 39-46, 8char
					k[6]= str('%8.3f' % (float(CG_coord[group][i][2])/10)).rjust(8)	#X 47-54, 8char
					fgro.write('%s%s%s%s%s%s%s\n' % (k[0],k[1],k[2],k[3],k[4],k[5],k[6]))
				#except:
				#	buffer = 0
			if int(l[8]) in terminal_residues:
				#writing terminal line
				fout.write("TER".ljust(6)+l[1]+l[2]+" "*(4)+l[4]+l[5]+l[6]+l[7]+l[8]+"\n")
		
		print ("Default box size is 10nm X 10nm X 10nm. Changes can be made directly in ", grofile," or using editconf.")
		fgro.write('%10.5f%10.5f%10.5f' % (10,10,10))
		fout.close(); fgro.close()
		return 1
	def nucAtomsType(self,pdbfile):
		#lading structure
		P = PDBParser(PERMISSIVE=0)
		structure = P.get_structure('test', pdbfile)

		#obtaining residues
		residues = structure.get_residues()
		r=[x.get_resname().lstrip().rstrip() for x in residues]

		#identifying unique nuclrotide residues
		temp = list()
		nuc_res = ["DA","DG","DC","DT","A","G","T","C","U"]
		for x in r:	
			if x in nuc_res and x not in temp:
				temp.append(x)
		r=temp

		#storing and counting deoxy nucleotides
		deoxy_count = [x for x in r if x.startswith("D")]
		#storing and counting oxy nucleotides
		oxy_count = [x for x in r if not x.startswith("D")]
		
		#chekcing total reasidues is euql to dexo and oxy ribose nucleotide
		assert(len(r)-len(deoxy_count)==len(oxy_count))

		#types of coarsed grain atoms
		a_type = ["B"+x[1] for x in deoxy_count]		#B as base represntator
		a_type = a_type + ["B"+x for x in oxy_count if "B"+x not in a_type] + ["XP"]

		#Adding sugar bead
		if len(deoxy_count)!=0:
			a_type = a_type + ["DS"+x[1] for x in deoxy_count]
		if len(oxy_count)!=0:
			a_type = a_type + ["RS"+x for x in oxy_count]
		na_type = len(a_type)
		return (na_type,a_type)
	def getContacts(self,pdbfile,nativefile,cutoff,atomname,strength, native_pairs):
		print (">>> Calculating intra-Nucelic-Acids contacts")

		cutoff=float(cutoff)
		#loading all atom pdb file
		load_atom = PDBParser().get_structure("test",pdbfile).get_atoms()
		#storing atom objects in a list
		atom = [a for a in load_atom]

		#identifying all atom contacts
		#Note: unlike protein CG models, nucleotide models should incorporate interactions between adj. residues (resi2-resi1>=1)

		#loading nativefile (CG model file)
		CG_file = PrePDB().readPDB(nativefile)
		CG_coord = dict()
		for l in CG_file:
			#key residue and atom number
			key = (int(l[22:26].strip()),l[12:16].strip())
			A = int(l[6:11])
			#X, Y, Z co-ordinates
			value = np.float_([l[30:38],l[38:46],l[46:54]])
			CG_coord[key] = [value,A]

		epsi_native = 1 #KJ/mol for native contacts
		CG_contacts = list()
		native_pairs = False #____TESTING_____
		if native_pairs:
			nuc_res = ("A","G","T","U","C","DA","DG","DT","DC")
			nuc_contacts = dict()
			for i in atom:
				for j in atom:
					#paret residue number
					r1 = i.get_parent().id[1]	
					r2 = j.get_parent().id[1]
					#parent residue name #only for re-checklng that the entries are intra nuc
					rname1 = i.get_parent().get_resname().strip()
					rname2 = j.get_parent().get_resname().strip()
					if r2!=r1 and rname1 in nuc_res and rname2 in nuc_res:
						d = float(j-i)		#distance
						if d <= cutoff:
							#atom number
							ai = i.get_id();aj = j.get_id()
							#residue number and atom number pair hashed to distance and residue anme
							if ((r2,r1),(aj,ai)) not in nuc_contacts:
								nuc_contacts[(r1,r2),(ai,aj)] = [d,(rname1,rname2)]
			print ("Total contacts within nucleiotides = ",len(nuc_contacts))

			#defining atoms in CG bead
			atoms_in = dict()
			atoms_in = PrePDB().atomIn()
			
			for pair,distance in nuc_contacts.items():
				pair_id  = pair[0]			#residue number pair
				pair_atom = pair[1]	#atomname pair
				pair_resname = distance[1]	#residue name

				a = dict()
				for ai in range(0,len(pair_atom)):
					#identifying Coarse grain atoms
					for group,atoms in atoms_in.items():
						if pair_atom[ai] in atoms:
							a[ai] = group
						else:
							continue	#skip this loop to next group as atom not in the group

					#assigning ribose type to the sugar bead
					if a[ai] == "S":
						if pair_resname[0][0] == "D":
							a[ai] = "D" + a[ai]	#adding D for deoxyribose
						else:
							a[ai] = "R" + a[ai]	#adding R for ribose
	
					#assigning base type to the base bead
					elif a[ai] == "B":
						if pair_resname[ai][0] == "D":
							a[ai] = a[ai] + pair_resname[ai][1]  #Adding A/T/G/C
						else:
							a[ai] = a[ai] + pair_resname[ai][0]  #Adding A/U/G/c
				
				#invoking CG coordinates from nativefile
				if (a[0][1]=="P" and a[1][1]=="S") or (a[0][1]=="S" and a[1][1]=="P"):
					if (int(pair_id[0])-int(pair_id[1]))**2 ==1:
						continue
						#Skipping the loop if P and S are adj,
				#storing CG coordinates from nativefile and calculating distances
				coord1 = CG_coord[pair_id[0],a[0]][0]
				coord2 = CG_coord[pair_id[1],a[1]][0]
				CG_dist = (np.sum((coord2-coord1)**2))**0.5 
				#storing atom id grom the nativefile
				atom1_id = CG_coord[pair_id[0],a[0]][1]
				atom2_id = CG_coord[pair_id[1],a[1]][1]
				
				#hcecking if pair not allready in the CG_contacts
				if (atom1_id,atom2_id,1,CG_dist, 1) not in CG_contacts:
					CG_contacts.append((atom1_id,atom2_id,1,CG_dist,epsi_native))	


		else:
			#used in experimental runs
			CG_coord = list(CG_coord.items()); CG_coord.sort()
			bases = list()
			for k,v in CG_coord:
				#K[0] = residue number #k[1] = atom name #v[0] = xyz #v[1] = atom number
				if k[1].startswith("B"):
					bases.append((k[1],v[1]))
			stack_dist = 3.6 #A 
			for i in range(0,len(bases)-1):
				#bases[i][1] = atom number while [0] = atom name
				CG_contacts.append((bases[i][1],bases[i+1][1],1, stack_dist, strength[bases[i][0]][bases[i+1][0]]))
		CG_contacts.sort()
		return CG_contacts
	def writeGroNonbondparams(self,topfilename,atoms_in_native):
		f = open(topfilename, "a")
		#f.write('%s\n' % ('[ nonbond_params ]'))
		#add non-bonded r6 term in sop-sc model for all non-bonded non-native interactions.
		#f.write('%s\n' % ('; i j  func sigma(c10)    eps(c12)'))

		temp_list = ["A","T","G","C","U","P","S"]
		basepairs= [("BA","BT"),("BA","DST"),("DSA","BT"),("DSA","DST"),
					("BA","BU"),("BA","RSU"),("RSA","BU"),("RSA","RSU"),
					("BG","BC"),("BG","DSC"),("DSG","BC"),("DSG","DSC"),
								("BG","RSC"),("RSG","BC"),("RSG","RSC")	]
		#_____TESTING______#
		basepairs=[]
		#_____TESTING______#
		

		for x in range(0,len(atoms_in_native)):
			i = atoms_in_native[x]
			for y in range(x,len(atoms_in_native)):
				j = atoms_in_native[y]
				#r  = ((i[5]/10 + j[5]/10)/2)**12
				r  = ((i[5]/10 + j[5]/10))**12 
				f.write('  %s\t%s\t1\t%e\t%e\n' % (i[0].strip(),j[0].strip(),0,r))
		f.close()
	def writeGroAtomTypes(self,topfilename,pdbfile,rad,no_Pcharge,excl_rule):
		print (">>> writing Nucleic Acids GROMACS topology atomtypes section",topfilename,pdbfile)
		#getting atom types
		(natomtypes,a_type) = self.nucAtomsType(pdbfile)

		#defining atomic unit mass for all Coarse grain atoms
		atommass=np.ones(natomtypes,dtype=np.float32)

		#defining charge
		atomcharge = dict()
		charge_in_atomtypes = False #to be modiified later
									#__TESTING__#
		for v in a_type:
			if not no_Pcharge and charge_in_atomtypes:
				if v == "XP":
					atomcharge[v] = -1
				else:
					atomcharge[v] = 0
			else:
				atomcharge[v] = 0
		assert(len(atomcharge)==natomtypes)
		
		#columns with known constant values
		c10 = 0.000
		atomptype = ["A"]*natomtypes

		#Assiging radius N-bases based on pyrimidine and purine input radius
		bases = dict()
		bases["Bpy"] = ("T","C","U")
		bases["Bpu"] = ("A","G")
		for group,base in bases.items():
			for b in base:
				rad[b] = rad[group]
		a = list()
		for i in range(0,len(a_type)):
			a.append([a_type[i],atommass[1],atomcharge[a_type[i]],atomptype[i],c10,rad[a_type[i][1]]])
		
		#opening topfile in append mode
		f = open(topfilename, "a")
		#writing section headrers
		f.write('%s\n'%('[ atomtypes ]'))
		f.write('%s\n' % ('; name mass  charge ptype c10    c12'))
		for i in (a):
			radius=(i[5]/10)**12 #C12 term #/10 for A to nm
			i[0]=i[0].rjust(5)
			f.write('%s %8.3f %8.3f %s\t%e\t%e \n'%(i[0],i[1],i[2],i[3],i[4],radius))


		#add non-bonded
		f.write('%s\n' % ('[ nonbond_params ]'))
		f.write('%s\n' % ('; i j  func sigma(c10)    eps(c12)'))
		f.close()
		if excl_rule == 2:
			#use arith. mean
			self.writeGroNonbondparams(topfilename,a)
		#self.writeBasePairs(topfilename,a)
		return a	
	def writeGroAtoms(self,topfilename, pdbfile, nativefile, no_Pcharge):
		print (">>> writing Nucleic Acids GROMACS topology atoms section")
		f1 = open(topfilename,"a")
		f1.write('\n%s\n'%('[ atoms ]'))
		f1.write('%s\n' % ("; nr  type  resnr  residue  atom  cgnr  charge  mass"))
		native_lines = PrePDB().readPDB(nativefile)
		l = {}		#atoms section line
		iter = 0	#atom number iterator
		for i in native_lines:
			iter = iter + 1
			l[0] = str(iter)				#atomnumber
			l[1] = i[12:16].strip()			#atom type
			if l[1].endswith("S"):			#adding residue type to the sugar
				l[1] = l[1] + i[17:20].strip()[len(i[17:20].strip())-1]
			l[2] = i[22:26].strip()			#residue number
			l[3] = i[17:20].strip()			#residue name
			l[4] = i[12:16].strip()			#atom name
			l[5] = str(iter)
			if i[12:16].lstrip().rstrip() == "XP":
				if not no_Pcharge:
					l[6] = str('%6.2f' % float(-1.000))
				elif no_Pcharge:
					l[6] = str('%6.2f' % float(0.000))
			else:
				l[6] = str('%6.2f' % float(0.000))
			l[7] = str('%6.2f' % float(1.000))
			f1.write("  ")
			for i in range(0,len(l)):
				f1.write(l[i]+"\t")
			f1.write("\n")
		f1.close()
		return 1
	def writeGroPairs(self,topfilename, nativefile, pdbfile, cutoff, atomname, native_pairs):
		print (">>> writing Nucleic Acids GROMACS pairs sections",topfilename, nativefile, pdbfile,cutoff)

		#Defining LJ-based stacking model Epsilon
		stackeps = PrePDB().stackEps()
		
		#get native conatcst #only if necessary
		CG_contacts=self.getContacts(pdbfile,nativefile,cutoff,atomname,stackeps,native_pairs)

		#writing in topfile
		f1 = open(topfilename, "a")
		f1.write('\n%s\n' % ('[ pairs ]'))
		f1.write('%s\n' % ('; ai\taj\ttype\tA\tB'))
		for i in CG_contacts:
			epsilonij=float(i[4])
			B=((i[3]*0.1)**12)*5*epsilonij	#C12 term
			A=((i[3]*0.1)**10)*6*epsilonij	#C10 term
			f1.write('  %d\t%d\t%s\t%e\t%e\n' % (i[0],i[1],'1',A,B))
		f1.close()
		
		bp_excl = self.write_dsDuplexBasePairs(topfilename,atomname,nativefile)
		f1 = open(topfilename, "a")
		f1.write('\n%s\n' % ('[ exclusions ]'))
		f1.write('%s\n' % ('; ai\taj'))

		#writing in exclusions
		if native_pairs:
			for i in CG_contacts:
				f1.write('  %d\t%d\n'%(i[0],i[1]))
		else:
			f1.write(";base_stacking\n")
			for i in CG_contacts:
				f1.write('  %d\t%d\n'%(i[0],i[1]))
		f1.write(";base_pairs;duplex\n")
		for i in bp_excl:
			f1.write('  %d\t%d\n'%(i[0],i[1]))
		f1.close()
		return 1
	def writeGroBonds(self,topfilename, nativefile, nKb, bond_ptype, CG_indices,terminal_residues):
		print (">>> writing Nucleic Acids GROMACS topology bonds section", topfilename, nativefile, nKb)
		#GROMACS 4.5.2 : FENE=7 AND HARMONIC=1
        #GROMACS IMPLEMENTS Ebonds = (Kx/2)*(r-r0)^2
        #Input units KJ mol-1 A-2
        #GROMACS units KJ mol-1 nm-1 (100 times the input value) 

		nKb = np.float_(nKb*100)	#KJ mol-1 nm-1


		bonds = list()
		for ai in range(1,len(CG_indices)):
			res1 = CG_indices[ai][0]		#residue number
			atm1 = CG_indices[ai][1]		#atom number
			xyz1 = CG_indices[ai][2]		#co-ordinates
			for aj in range(ai+1,len(CG_indices)+1):
				res2 = CG_indices[aj][0]	#residue numbeer
				atm2 = CG_indices[aj][1]	#atom number
				xyz2 = CG_indices[aj][2]	#co-ordintaes
				if res1 in terminal_residues and 0 > res1-res2:
					#skip loop if res1 is terminal residue and
					#res2 is from the next chain
					continue
				if res2-res1 == 0: 
					#Bonding term for same residue 
					#Pi-Si
					if (atm1[1]=="P" and atm2[1]=="S") or (atm1[1]=="S" and atm2[1]=="P"):
						bonds.append((ai,aj,(np.sum((xyz2-xyz1)**2))**0.5))
					#Si-Bi
					elif (atm1[1]=="S" and atm2[0]=="B") or (atm1[0]=="B" and atm2[1]=="S"):
						bonds.append((ai,aj,(np.sum((xyz2-xyz1)**2))**0.5))
				elif res2-res1 == 1:
					#Bonding term for adj resi 
					#Si-Pi+1
					if (atm1[1]=="S" and atm2[1]=="P"):
						bonds.append((ai,aj,(np.sum((xyz2-xyz1)**2))**0.5))
		#writing topology bonds section
		f1 = open(topfilename, "a")
		f1.write('\n%s\n' % ('[ bonds ]'))
		f1.write('%s\t%s\t%s\t%s\t%s\n' % ('; ai', 'aj', 'func', 'r0(nm)', 'Kb'))
		#force constant
		bonds.sort()
		for i in bonds:
			f1.write('  %d\t%d\t%d\t%e\t%e\n' % (i[0], i[1], bond_ptype, i[2]/10, nKb))
		f1.close()
		return bonds
	def writeGroAngles(self,topfilename, nativefile, nKa, CG_indices,terminal_residues):
		print (">>> writing Nucleic Acids GROMACS topoplogy angle section", topfilename, nativefile, nKa)
        #Eangles = (Ktheta/2)*(r-r0)^2
        #Input units KJ mol-1
        #GROMACS units KJ mol-1  
		nKa = float(nKa)

		triplets = list()
		#segregating indices based on atom
		resi_dict = dict()
		for key,value in CG_indices.items():
			resi = value[0]	#residue number
			atom = value[1]	#atom name
			xyz  = value[2]	#co-ordinates
			if resi not in resi_dict:
				resi_dict[resi] = dict()
			if atom[1] in ("P","S"):
				resi_dict[resi][atom[1]] = (key,xyz)
			elif atom[0] == "B":
				resi_dict[resi][atom[0]] = (key,xyz)

		triplets = list()
		for r in range(1,len(resi_dict)):
			#r here is the residue number
			try:#Set I P(i)S(i)B(i)
				triplets.append((resi_dict[r]["P"][0]-1,resi_dict[r]["S"][0]-1,resi_dict[r]["B"][0]-1))	#making indices start from 0
			except:
				buf = 0
			if r not in terminal_residues:
				#if i is terminal residue, don't include the triplet
				try:#Set II B(i)S(i)P(i+1)
						triplets.append((resi_dict[r]["B"][0]-1,resi_dict[r]["S"][0]-1,resi_dict[r+1]["P"][0]-1))	#making indices start from 0
				except:
					buf = 0
				try:#Set III P(i)S(i)P(i+1)
					triplets.append((resi_dict[r]["P"][0]-1,resi_dict[r]["S"][0]-1,resi_dict[r+1]["P"][0]-1))	#making indices start from 0
				except:
					buf = 0
				try:#Set IV S(i)P(i+1)S(i+1)
					triplets.append((resi_dict[r]["S"][0]-1,resi_dict[r+1]["P"][0]-1,resi_dict[r+1]["S"][0]-1))	#making indices start from 0
				except:
					buf=0

		M = Measure()
		CG_angles = M.angles(nativefile,triplets)
		CG_angles = CG_angles[0]

		angles = list()
		for i in range(0,len(CG_angles)):
			 angles.append(((triplets[i]+np.array([1,1,1])).tolist(),CG_angles[i]))
		angles.sort()

		f1 = open(topfilename, "a")
		f1.write('\n%s\n' % ('[ angles ]'))
		f1.write('%s\t%s\t%s\t%s\t%s\t%s\n' % ('; ai', 'aj', 'ak','func', 'th0(deg)', 'Ka'))
		for i in angles:
			f1.write('  %d\t%d\t%d\t%s\t%e\t%e\n' % (i[0][0],i[0][1],i[0][2],'1',i[1]*180/np.pi, nKa))
		f1.close()
		return angles
	def writeGroDihedrals(self,topfilename, nativefile, base_Kd, backbone_Kd, CG_indices,terminal_residues,P_stretch):
		print (">> writing GROMACS dihedrals section", topfilename,nativefile,"Base Kd",base_Kd,"backbone Kd",backbone_Kd)
		#converting Kcal/mol to KJ/mol
		#K_dihedrals will be divided by these factprs
		#_____USER_CAN_EDIT____#
		#The Kd will be multiplied the factor below
		factor_1 = 1.0		#for term with multiplicity 1
		factor_2 = 1.0		#for term with multiplicity 3
		#______________________#

		#GROMACS IMPLEMENTATION: Edihedrals Kphi*(1 + cos(n(phi-phi0)))
		#Our implementaion: Edihedrals = Kphi*(1 - cos(n(phi-phi0)))
		#The negative sign is included by add phase = 180 to the phi0
		#Kphi*(1 + cos(n(phi-180-phi0))) = Kphi*(1 + cos(n180)*cos(n(phi-phi0)))
		#if n is odd i.e. n=1,3.... then cos(n180) = -1
		#hence Edihedrals = Kphi*(1 - cos(n(phi-phi0)))

		phase = 180

		
		
		quadruplets_base = list()
		quadruplets_backbone = list()
		resi_dict = dict()
		for key,value in CG_indices.items():
			resi = value[0]
			atom = value[1]
			xyz  = value[2]
			if resi not in resi_dict:
				resi_dict[resi] = dict()
			if atom[1] in ("P","S"):
				resi_dict[resi][atom[1]] = (key,xyz)
			elif atom[0] == "B":
				resi_dict[resi][atom[0]] = (key,xyz)
		
		for r in range(1,len(resi_dict)):
			if r not in terminal_residues:
				try: #Set I B(i)S(i)S(i+1)B(i+1)
					quadruplets_base.append((resi_dict[r]["B"][0]-1,resi_dict[r]["S"][0]-1,resi_dict[r+1]["S"][0]-1,resi_dict[r+1]["B"][0]-1))
				except:
					buf = 0		#buffer element for exception
			if r not  in terminal_residues and r+1 not in terminal_residues and r+2 not in terminal_residues:
				#print (r,r+1,r+2,r+3)
				try: #Set II P(i)P(i+1)P(i+2)P(i+3)				
					quadruplets_backbone.append((resi_dict[r]["P"][0]-1,resi_dict[r+1]["P"][0]-1,resi_dict[r+2]["P"][0]-1,resi_dict[r+3]["P"][0]-1))
				except:
					buf = 0
		
		M = Measure()
		CG_base_dihedrals = M.dihedrals(nativefile,quadruplets_base)
		CG_base_dihedrals = CG_base_dihedrals[0]
		
		#backbone dihedrals are extended to 180 degreees so no need to aquire
		#CG_backbone_dihedrals = Y.get_dihedrals(nativefile,quadruplets_backbone)[0]
		#___TESTING____
		CG_backbone_dihedrals = M.dihedrals(nativefile,quadruplets_backbone)[0]

		f1 = open(topfilename, "a")
		f1.write('\n%s\n' % ('[ dihedrals ]'))
		f1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('; ai','aj','ak','al','func','phi0(deg)','Kd','mult'))

		dihedrals = list()
		for i in range(0,len(quadruplets_base)):
			dihedrals.append(((quadruplets_base[i]+np.array([1,1,1,1])).tolist(),CG_base_dihedrals[i], base_Kd))
		for i in range(0,len(quadruplets_backbone)):
			if P_stretch:
				dihedrals.append(((quadruplets_backbone[i]+np.array([1,1,1,1])).tolist(),np.pi, backbone_Kd))
			else:
				dihedrals.append(((quadruplets_backbone[i]+np.array([1,1,1,1])).tolist(),CG_backbone_dihedrals[i], backbone_Kd))
		dihedrals.sort()

		for i in dihedrals:
			#dividing by factor_1 and factor_2 for 1 and 3 respectively
			#multiplicity 1
			phi0 = phase+i[1]*180/np.pi	#adding 180 degree to hcange phase
			f1.write(' %d\t%d\t%d\t%d\t%s\t%e\t%e\t%s\n'% (i[0][0], i[0][1], i[0][2], i[0][3], '1 ',1*phi0, i[2]/factor_1, "1"))
			#multiplicity 3				#using 3*phi0
			f1.write(' %d\t%d\t%d\t%d\t%s\t%e\t%e\t%s\n'% (i[0][0], i[0][1], i[0][2], i[0][3], '1 ',3*phi0, i[2]/factor_2, "3"))
		f1.close()
		return dihedrals
	def writeBasePairs(self,topfilename,atoms_in_native):
		#print (">> writing GROMACS pairtypes section section", topfilename, nativefile, nKa)
		
		atoms_in_native = [x[0].strip() for x in atoms_in_native]	#x[0] is atom types
		ribose = 0; deoxy_ribose = 0
		for i in atoms_in_native:
			if i.startswith("D"):
				print (i)
				deoxy_ribose = 1
			elif i.startswith("R"):
				ribose = 1

		if ribose == 1 and deoxy_ribose == 0:
			inputfile = "a.rna.pdb"
		elif deoxy_ribose == 1 and ribose == 0:
			inputfile = "b.dna.pdb"
		elif ribose + deoxy_ribose == 2:
			inputfile = "h_rna_dna.pdb"
		sample_file = open(inputfile)

		#loading bead co-orifnates for the sample file
		#The sample file is generated using the same paramaters as defined by the 
		#user for the model system. PDB co-ordinates are present in the function itself

		sample = [x for x in sample_file if x.startswith("ATOM")]
		bases = [x for x in sample if x[12:16].strip().startswith("B")]
		sugars = [x for x in sample if x[12:16].strip().endswith("S")]

		#defining LJ strenfth for base pairs (will be effectively 4 times of the below value)
		epsilon = dict()
		epsilon = PrePDB().basepairEps()

		#loading WC-pairing dist from the sample file
		h_bond_dist = dict()
		backbone_dist =  dict()
		sug_base_dist_1 = dict()
		sug_base_dist_2 = dict()
		rescount = len(bases)
		for i in range(0,rescount/2):
			#loading pdb lines
			base1 = bases[i].strip()
			base2 = bases[rescount-1-i].strip()
			sugar1 = sugars[i].strip()
			sugar2 = sugars[rescount-1-i].strip()
						
			#loading atom names
			base_atom1 = base1[12:16].strip()
			base_atom2 = base2[12:16].strip()

			#sugar atom type (ending with the single letter code for the nucleotide)
			sugar_atom1 = sugar1[12:16].strip() + base_atom1[1]	#for BN atom sugar type is RSN or DSN
			sugar_atom2 = sugar2[12:16].strip() + base_atom2[1] #for BN atom sugar type is RSN or DSN

			#loading co-ordinates
			#bases
			xyz_b1 = np.float_([base1[30:38],base1[38:46],base1[46:54]])
			xyz_b2 = np.float_([base2[30:38],base2[38:46],base2[46:54]])
			#sugars
			xyz_s1 = np.float_([sugar1[30:38],sugar1[38:46],sugar1[46:54]])
			xyz_s2 = np.float_([sugar2[30:38],sugar2[38:46],sugar2[46:54]])
			
			#H-bond base-base distance
			h_bond_dist[i] = [0.1*(sum((xyz_b2-xyz_b1)**2))**0.5,base_atom1,base_atom2] #0.1 multipled to convert A to nm
			#Base pairing distance b/w sugars of pairing bases

			backbone_dist[i] = [0.1*(sum((xyz_s2-xyz_s1)**2))**0.5,sugar_atom1,sugar_atom2] #0.1 multiples to convert A to nm
			
			#base-sugar and sugar-base terms generic contacts
			sug_base_dist_1[i] = [0.1*(sum((xyz_b2-xyz_s1)**2)**0.5),sugar_atom1,base_atom2]
			sug_base_dist_2[i] = [0.1*(sum((xyz_b1-xyz_s2)**2)**0.5),base_atom1,sugar_atom2]

		fout = open(topfilename, "a")
		#header = "nonbond_params"
		#fout.write("\n[ "+header+" ]\n")
		#fout.write("; ai\taj\tfunc\tC10\tC12\n")
		fout.write(";base_pairs" + "\n")

		#pre-defined list of to-be-added parameteres
		parameters = [h_bond_dist,sug_base_dist_1,sug_base_dist_2,backbone_dist]
		ftest = open("base_paris.t","w+")
		
		for distance  in parameters:
			#write base1-base2
			#wrtie sugar1-base2 and sugar2-base1
			#write sugar1-sugar2 distacne
			for k in range(len(distance)):
				#atom name stored at inedx 0 and 1
				atom1 = distance[k][1]
				atom2 = distance[k][2]
				if atom1 not in atoms_in_native or atom2 not in atoms_in_native:
					#skip, the pair doesnot exist in the structure itself
					continue	#moving to next pair

				#value stored at index 00
				value = distance[k][0]

				#to check if the atom exists in the native file
				#if  atom1 not in atoms_in_native or atom2 not in atoms_in_native:
				#	continue
				
				#for all above atom type, the string ends with nucleotide 1 letter code
				strength = epsilon[(atom1[len(atom1)-1],atom2[len(atom2)-1])]
				
				B=(value**12)*5*strength
				A=(value**10)*6*strength
				
				print (atom1,"\t",atom2,"\t",value,"\n")
				ftest.write('%s\t%s\t%e\n' % (atom1,atom2,value))
				fout.write('  %s\t%s\t%d\t%e\t%e\n' % (atom1,atom2,1,A,B))
		ftest.close()
		fout.close()
		return 1
	def write_dsDuplexBasePairs(self,topfilename,atoms_in_native,nativefile):
		
		#The function generates all possible base pairing duplex pairs. 
		# Total unique ds-duplex = 10 (exclusing wobble G-U pair) 
		# The code will check for prsence of each unique ds-duplex in the sequence
		# and a pair will be defined wbetween corrosponding nucleotides for 
		# duplex pair match. 
		# Eg-> A and G 5'-AG-3' will paur wilh U and C of 5'-CU-3' respectively
		  	
		atoms_in_native = [x[0].strip() for x in atoms_in_native]	#x[0] is atom types
		ribose = 0; deoxy_ribose = 0
		for i in atoms_in_native:
			if i.startswith("D"):
				print (i)
				deoxy_ribose = 1
			elif i.startswith("R"):
				ribose = 1

		if ribose == 1 and deoxy_ribose == 0:
			inputfile = "a.rna.pdb"
			#complementry nucleotide
			c = {'A':'U','G':'C','U':'A','C':'G'}
		elif deoxy_ribose == 1 and ribose == 0:
			inputfile = "b.dna.pdb"
			#complementry nucleotide
			c = {'A':'T','G':'C','T':'A','C':'G'}
		elif ribose + deoxy_ribose == 2:
			inputfile = "h_rna_dna.pdb"
		sample_file = open(inputfile)

		#loading bead co-orifnates for the sample file
		#The sample file is generated using the same paramaters as defined by the 
		#user for the model system. PDB co-ordinates are present in the function itself

		sample = [x for x in sample_file if x.startswith("ATOM")]
		bases = [x for x in sample if x[12:16].strip().startswith("B")]
		sugars = [x for x in sample if x[12:16].strip().endswith("S")]

		#defining LJ strenfth for base pairs (will be effectively 4 times of the below value)
		epsilon = dict()
		epsilon = PrePDB().basepairEps()

		#loading WC-pairing dist from the sample file
		h_bond_dist = dict()
		backbone_dist =  dict()
		sug_base_dist = dict()
		rescount = len(bases)
		for i in range(0,int(rescount/2)):
			#loading pdb lines
			base1 = bases[i].strip()
			base2 = bases[rescount-1-i].strip()
			sugar1 = sugars[i].strip()
			sugar2 = sugars[rescount-1-i].strip()
						
			#loading atom names
			base_atom1 = base1[12:16].strip()
			base_atom2 = base2[12:16].strip()

			#sugar atom type (ending with the single letter code for the nucleotide)
			sugar_atom1 = sugar1[12:16].strip() + base_atom1[1]	#for BN atom sugar type is RSN or DSN
			sugar_atom2 = sugar2[12:16].strip() + base_atom2[1] #for BN atom sugar type is RSN or DSN

			#loading co-ordinates
			#bases
			xyz_b1 = np.float_([base1[30:38],base1[38:46],base1[46:54]])
			xyz_b2 = np.float_([base2[30:38],base2[38:46],base2[46:54]])
			#sugars
			xyz_s1 = np.float_([sugar1[30:38],sugar1[38:46],sugar1[46:54]])
			xyz_s2 = np.float_([sugar2[30:38],sugar2[38:46],sugar2[46:54]])
			
			#H-bond base-base distance
			h_bond_dist[(base_atom1,base_atom2)] = 0.1*(sum((xyz_b2-xyz_b1)**2))**0.5 #0.1 multipled to convert A to nm
			h_bond_dist[(base_atom2,base_atom1)] = h_bond_dist[(base_atom1,base_atom2)] 
			#Base pairing distance b/w sugars of pairing bases

			backbone_dist[(sugar_atom1,sugar_atom2)] = 0.1*(sum((xyz_s2-xyz_s1)**2))**0.5 #0.1 multiples to convert A to nm
			backbone_dist[(sugar_atom2,sugar_atom1)] = backbone_dist[(sugar_atom1,sugar_atom2)] 

			#base-sugar and sugar-base terms generic contacts
			sug_base_dist[(sugar_atom1,base_atom2)] = 0.1*(sum((xyz_b2-xyz_s1)**2)**0.5)
			sug_base_dist[(base_atom2,sugar_atom1)] = sug_base_dist[(sugar_atom1,base_atom2)]
			sug_base_dist[(sugar_atom2,base_atom1)] = 0.1*(sum((xyz_b1-xyz_s2)**2)**0.5)
			sug_base_dist[(base_atom1,sugar_atom2)] = sug_base_dist[(sugar_atom2,base_atom1)]

		#reading natiev file to get pdb sequence and co-ordinates
		fin = open(nativefile)
		seq = ""
		base = {}
		sugar = {}
		phos = {}
		base_name = {}
		sugar_name = {}
		for i in fin:
			if i.startswith("ATOM"):
				resnum = int(i[22:26])
				atnum = int(i[6:11])
				resname = i[17:20].strip()[len(i[17:20].strip())-1]
				atname = i[12:16].strip()
				if atname[0] == "B":
					#reading sequence
					seq=seq+atname[1]
					#adding atunm number
					base[resnum] = atnum
					base_name[resnum] = atname
				elif "P" in atname:
					phos[resnum] = atnum
				elif "S" in atname:
					sugar[resnum] = atnum
					sugar_name[resnum] = atname + resname
		#list of nucleotide in the 
		n = list()
		for i in seq:
			if i not in n:
				n.append(i)
		fin.close()



		#Generating unique duplex pairs#
		duplex = list()
		for i in n:
			for j in n:
				if (i+j,c[j]+c[i]) not in duplex and (c[j]+c[i],i+j) not in duplex:
					#pair if not in duplex list.
					#checking both sides in 5'->3' order
					duplex.append((i+j,c[j]+c[i]))
		print ("number of Unique pairs: ",len(duplex))
		print ("unique pairs: ",duplex)
		
		#searching position of each index in the sequence
		position = dict()
		for pair in duplex:
			for p in pair:
				if p not in position:
					#creating a position list for duplex p
					position[p] = list()
				index = 0 #initiating index
				for i in seq.split(p):
					#spliting sequence at residue p and adding the fregment length
					index = index + len(i)
					if index >= len(seq):
						continue	#skip if index exceeds the seq length
					#appending index to the list
					position[p].append(index+1)
					index = index + len(p)	#adding length of duplex

		#removing ds-duplex with not occurance in the seq
		temp = list()
		for pair in duplex:
			print (seq,pair,len(position[pair[0]]),len(position[pair[1]]))
			if len(position[pair[0]]) == 0 or len(position[pair[1]]) == 0:
				continue
			#add pair to temp list if the ds-duplex exists in the seq
			temp.append(pair)
		duplex=temp	#replace input by filtered output
		del(temp)

		#generating pair list
		temp = dict() #only for checks and will be deleted
		plist = dict()
		for p,q in duplex:
			temp[(p,q)] = list()
			plist[(p,q)] = list()
			for i in position[p]:		#index of p
				for j in position[q]:	#index of q
					if i!=j and ((i-j)**2)**0.5 > 2:
					#not pairing to iteself or to adj. residue
						if (i,j) not in temp[(p,q)] and (j,i) not in temp[(p,q)]:
							#checking if pair not allready defined
							if len(p)+len(q) == 4:
								#for ds-duplex
								temp[(p,q)].append((i,j))
								# AG-CU (5AG3-3UC5) duplex appending A-u pair and GC pair
								plist[(p,q)].append(((i,j+1),(i+1,j)))
							if len(p)+len(q) == 6:
								temp[(p,q)].append((i,j))
								plist[(p,q)].append(((i,j+2),(i+2,j),(i+1,j+1)))
		duplex = plist
		del(temp); del(plist)

		#writing pairs in pairs section
		fout = open(topfilename,"a")
		fout.write(";;;;;base-pairing;;;ds-duplex;;;;;")
		excl = list()

		for basepair in duplex:
			for set in duplex[basepair]:
				for i in set:
					D = h_bond_dist[(base_name[i[0]],base_name[i[1]])]
					excl.append((base[i[0]],base[i[1]]))
					fout.write ('\n  %d\t%d\t%d\t%e\t%e\t' % (base[i[0]],base[i[1]],1,6*D**10,5*D**12))
					D = sug_base_dist[(base_name[i[0]],sugar_name[i[1]])]
					excl.append((base[i[0]],sugar[i[1]]))
					fout.write ('\n  %d\t%d\t%d\t%e\t%e\t' % (base[i[0]],sugar[i[1]],1,6*D**10,5*D**12))
					D = sug_base_dist[(sugar_name[i[0]],base_name[i[1]])]	
					excl.append((sugar[i[0]],base[i[1]]))
					fout.write ('\n  %d\t%d\t%d\t%e\t%e\t' % (sugar[i[0]],base[i[1]],1,6*D**10,5*D**12))
					D = backbone_dist[(sugar_name[i[0]],sugar_name[i[1]])]
					excl.append((sugar[i[0]],sugar[i[1]]))
					fout.write ('\n  %d\t%d\t%d\t%e\t%e\t' % (sugar[i[0]],sugar[i[1]],1,6*D**10,5*D**12))
		fout.close()
		return excl
	def writeGromacsTop(self,topfilename,pdbfile,nativefile,btparams,rad,fconst,cutoff,no_Pcharge,native_pairs,P_stretch,excl_rule):
		#cheking if nucleotide exists
		checksize = open(pdbfile)
		if len(checksize.readlines()) <= 1:
			#if pdbfile size = 0, then there is no nucleotides im the sysytem
			#exiting with 0
			#print ("As mentioned earlier, the ibpu file is PROTEIN ONLY.....Exiting now")
			return 0
		checksize.close()

		print ("====>>>>> Generating Nucleotide CG files <<<<=========")
		#Order for SBM file
		#2 = 10-12 potential
		# Note: Parameters in GROMACS units need tweaking.
		#write gromacs format files.
		CG_indices,terminal_residues = self.getIndices(nativefile)


		bond_ptype = 1
		self.write_gro_header(topfilename,3) 	#3 is atom type 		#inheriting from protSBM
		self.write_header_SBM()						  		#inheriting from protSBM	
		atomname = self.writeGroAtomTypes(topfilename,pdbfile,rad,no_Pcharge,excl_rule)
		self.write_gro_moleculetype(topfilename)
		self.writeGroAtoms(topfilename, pdbfile, nativefile, no_Pcharge)
		self.writeGroBonds(topfilename, nativefile, fconst["nKb"], bond_ptype, CG_indices,terminal_residues)
		self.writeGroAngles(topfilename, nativefile, fconst["nKa"], CG_indices,terminal_residues)
		self.writeGroDihedrals(topfilename, nativefile, fconst["nKd"], fconst["P_nKd"], CG_indices,terminal_residues,P_stretch)
		self.writeGroPairs(topfilename, nativefile, pdbfile, cutoff, atomname, native_pairs)
		self.write_gro_tail(topfilename)
		print ('See file:',topfilename,'for GROMACS topology.')
		return 1
	def writeTablefile(self,debye,D,iconc,irad,T,power):
		print (">>Print table file",'table_debye_'+str(power)+'_12.xvg')
		pi = np.pi
		irad = irad/10		#convertig A to nm
		jtoc = 0.239		#C/J
		if debye:
			#debye length
			#	 			e0*D*KB*T
			#dl**2 = -------------------------
			#		 2*el*el*NA*l_to_nm3*iconc


			e0 = 8.854e-21 			#C2 J-1 nm-1	#permitivity: The value have been converted into C2 J-1 nm-1
			D = D									#dielectric constant D
			KB = 1.3807*10**(-23) 	#J/K			#Boltzmann constant
			el = 1.6e-19			#C				#charge on electron
			NA = 6.022e+23          #n/mol  		#Avogadro's number
			T = 298					#K				#Temperature
			l_to_nm3 = 1e-24		#nm3/l			#converting l to nm3

			dl_N = e0*D*KB*T	#numerator
			dl_D = 2*el*el*NA*l_to_nm3		#denom
			dl = (dl_N/dl_D)**0.5 #debye length nm M L-1
			#Note: still testing
			iconc = iconc 
			dl = 1/dl			#nm-1 M-1 L
			dl = dl*(iconc**0.5) 		#nm-1
			Bk = np.exp(dl*irad)/(1+dl*irad)
		
		else:
			dl = 0
			Bk = 1
		Ke = 1/(4*pi*8.854e-21)                     #J nm C-2
		Ke = Ke*(2.56e-38)*(6.022e+23)*(1e-3)       #KJ nm C-2 e-2
		Ke=jtoc*Ke*Bk/(D)
		#Joule_Cal_jhol
		
		fout = open('table_debye_'+str(power)+'_12.xvg','w+')
		fout1 = open("table_nuc-aa_RNA_RNA.xvg","w+")
		for i in range(0,2*29018+2):
			r = float(i)*0.002
			if r > 0.01:
				if r>=0.01:
					V 		= jtoc*(Bk/D)*np.exp(-dl*r)/r
					V_1 	= jtoc*(Bk/D)*(-np.exp(-dl*r)/r**2) + (Bk/D)*(-dl*np.exp(-dl*r)/r)
					V_1		= V_1*(-1)
				else:
					V = 0
					V_1 = 0
				Cn 	= -1/r**power 	
				Cn_1	= -10/r**(power+1) 
				C12		= 1/r**12
				C12_1	= 12/r**13
				fout.write('%e %e %e %e %e %e %e\n' %(r,V,V_1,Cn,Cn_1,C12,C12_1))
				fout1.write('%e %e %e %e %e %e %e\n' %(r,0,0,Cn,Cn_1,C12,C12_1))
			else:
				fout.write('%e %e %e %e %e %e %e\n' %(r,0,0,0,0,0,0))
				fout1.write('%e %e %e %e %e %e %e\n' %(r,0,0,0,0,0,0))
		fout.close()
		fout1.close()
		return 	'table_debye_'+str(power)+'_12.xvg'

def main():
	import argparse
	parser = argparse.ArgumentParser(description="Generate GROMACS and OPTIM potential files for enhanced SBM models.")
	parser.add_argument("--P_rad", help="Radius for Backbone Phosphate group bead. Default=3.7A")
	parser.add_argument("--S_rad", help="Radius for Backbone Sugar group bead. Default=3.7A")
	parser.add_argument("--Bpu_rad", help="Radius for N-Base Purine bead. Default=1.5A")
	parser.add_argument("--Bpy_rad", help="Radius for N-Base Pyrimidine bead. Default=1.5A")
	parser.add_argument("--nKb", help="Kbond for RNA/DNA")
	parser.add_argument("--nKa", help="Kangle for RNA/DNA. Default=20")
	parser.add_argument("--nKd", help="Kdihedral for Bi-Si-Si+1-Bi+1. Default=0.5")
	parser.add_argument("--P_nKd", help="Kdihedral for Backbone Pi-Pi+1-Pi+2-Pi+3. Default=0.7")
	parser.add_argument("--cutoff", help="Cut-off for contact-map generation")
	
	parser.add_argument("--Bpu_pos", help="Put input atom of Purine [N1,C2,H2-N2,N3,C4,C5,C6,O6-N6,N7,C8,N9,COM] as position of B. Default=COM(Center_of_Mass)")
	parser.add_argument("--Bpy_pos", help="Put input atom of Pyrimidine [N1,C2,O2,N3,C4,O4-N4,C5,C6,COM] as position of B. Default=COM(Center_of_Mass)")
	parser.add_argument("--S_pos", help="Put input atom of Sugar [C1',C2',C3',C4',C5',H2'-O2',O3',O4',O5',COM] as position of S   . Default=COM(Center_of_Mass)")
	parser.add_argument("--P_pos", help="Put input atom of Phosphate [P,OP1,OP2,O5',COM] group as position of P. Default=COM(Center_of_Mass)")
	
	parser.add_argument("--scaling", help="Scaling for mapping to all-atom contact-map.")
	parser.add_argument("--btmap", action='store_true', default=False, help='Use Betancourt-Thrumalai parameters for Bead interactions.')
	parser.add_argument("--w_native", help='Write native files, P-S-B from all atom PDB file.')
	parser.add_argument("--grotop",help='Gromacs topology file output name.')
	parser.add_argument("--pl_map", action='store_true', default=False, help='Plot contact map for two bead model')
	parser.add_argument("--aa_pdb", help='all-atom pdbfile e.g. 1qys.pdb')
	parser.add_argument("--all_chains", action='store_true', default=False, help='Will Not remove identical chains.')	
	parser.add_argument("--pdb2gro", help='Name for .gro file.')
	parser.add_argument("--pistacklen", help="pi-pi stacking length. Default=3.6A")

	parser.add_argument("--debye",action='store_true', help="Use debye electrostatic term")
	parser.add_argument("--T", help="System temperature. Default = 298K")
	parser.add_argument("--Kr", help="Krepulsion. Default=5.7A")
	parser.add_argument("--iconc", help="Solvant ion conc. Default=0.1M")  
	parser.add_argument("--irad", help="Solvant ion rad. Default=1.4A")  
	parser.add_argument("--dielec", help="Dielectric constant of solvant. Default=70")
	parser.add_argument("--no_Pcharge","-no_Pcharge", action='store_true', default=False, help='No negative charge on Phosphate bead')
	parser.add_argument("--P_stretch",help="Stretch the backbone dihedral to 180 degrees. Default = Use native  backbone dihedral")

	#parser.add_argument("--base_pair","-base_pair", action='store_true', default=False, help='Use basepairing terms')
	parser.add_argument("--excl_rule",help="Use 1: Geometric mean. 2: Arithmatic mean")
	parser.add_argument("--Kr", help="Krepulsion. Default=5.7A")


	args = parser.parse_args()
	N = nucsbm()
	P_stretch = False
	#Set default parameters
	native_pairs = True
	no_Pcharge = False
	fconst = dict()					#force_csontant
	rad = dict()					#CG bead radii
	fconst["nKb"]= 200				#KCal/mol
	fconst["nKa"] = 20				#KCal/mol
	fconst["nKd"] = 0.5				#KCal/mol
	fconst["P_nKd"] = 0.7			#KCal/mol
	fconst["Kr"] = 5.7				#KCal/mol
	rad["P"] = 3.7					#A
	rad["S"] = 3.7					#A
	rad["Bpy"] = 1.5				#A
	rad["Bpu"] = 1.5				#A
	rad["stack"] = 3.6					#A
	cutoff = 4.5					#A
	irad = 1.4						#A (Na+ = 1.1A	Cl- = 1.7A)
	iconc = 0.01					#M
	D =	70							#dilectric constant (dimensionlsess)
	debye = False
	T = 120 #K
	pur_atom = ("N1","C2","H2-N2","N3","C4","C5","C6","O6-N6","N7","C8","N9","COM")
	pyr_atom = ("N1","C2","O2","N3","C4","O4-N4","C5","C6","H7-C7","COM")
	sug_atom = ("C1'","C2'","C3'","C4'","C5'","H2'-O2'","O3'","O4'","O5'","COM")
	phos_atom = ("P","OP1","OP2","O5'","COM")
	Bpu_pos = "COM"		#Center of Mass for purine
	Bpy_pos = "COM"		#Center of Mass for pyrimidine
	P_pos = "COM"			#Center of Mass for phosphate group
	S_pos = "COM"			#Center of Mass for sugar
	btparams = False	

	#Replacing default paramteres with input paramters
	if args.P_stretch:
		P_stretch = True
	if args.no_Pcharge:
		no_Pcharge = True
	if args.Bpu_pos:
		argsbuf = str(args.Bpu_pos)
		if argsbuf in pur_atom:
			Bpu_pos = argsbuf
		else:
			print("Warning!!! Wrong Atom name entered for Purine. The program will continue with default parameter")
	if args.Bpy_pos:
		argsbuf = str(args.Bpy_pos)
		if argsbuf in pyr_atom:
			Bpy_pos = argsbuf
		else:
			print("Warning!!! Wrong Atom name entered for Pyrimidine. The program will continue with default parameter")
	if args.P_pos:
		argsbuf = str(args.P_pos)
		if argsbuf in phos_atom:
			P_pos = argsbuf
		else:
			print("Warning!!! Wrong Atom name entered for Phosphate. The program will continue with default parameter")
	if args.S_pos:
		argsbuf = str(args.S_pos)
		if argsbuf in sug_atom:
			S_pos = argsbuf
		else:
			print("Warning!!! Wrong Atom name entered for Sugar. The program will continue with default parameter")
	if args.excl_rule:
		excl_rule = int(args.excl_rule)
		if excl_rule not in (1,2):
			print ("Choose correct exclusion rule. Use 1: Geometric mean or 2: Arithmatic mean")
	else:
		excl_rule = 1


	CG_pos = {"Bpu":Bpu_pos,"Bpy":Bpy_pos,"S":S_pos,"P":P_pos}	
	print (CG_pos)
	#Force constants
	if args.nKb:
	    fconst["nKb"] = float(args.nKb)
	if args.nKa:
		fconst["nKa"] = float(args.nKa)
	if args.nKd:
		fconst["nKd"] = float(args.nKd)
	if args.P_nKd:
		fconst["P_nKd"] = float(args.P_nKd)

	if args.Kr:
		fconst["Kr"] = float(args.Kr)
	#CG radii
	if args.P_rad:
		rad["P"] = float(args.P_rad)
	if args.S_rad:
		rad["S"] = float(args.S_rad)
	if args.Bpu_rad:
		rad["Bpu"] = float(args.Bpu_rad)
	if args.Bpy_rad:
		rad["Bpy"] = float(args.Bpy_rad)
	if args.pistacklen:
		rad["stack"] = float(args.pistacklen)
	
	#solvant and ionic params
	if args.dielec:
		D = float(args.dielec)
	if args.iconc:
		iconc = float(args.iconc)
	if args.irad:
		irad = float(args.irad)
	sol_params = {"conc":iconc,"rad":irad,"D":D}
	if args.debye:
		debye = True
	if args.T:
		T = float(args.T)
	if args.cutoff:
	    cutoff=float(args.cutoff)
	if args.btmap:
		btparams=True
	
	#cheking for read-write/input-output options
	grofile = ""
	if args.pdb2gro:
		grofile = str(args.pdb2gro)
	if args.w_native:
		pdbfile = args.w_native
		if not args.all_chains:
			pdbfile = PrePDB().homoNmer2Monomer(pdbfile)
		(pdbfile,terminal_residues) = PrePDB().appendChain(pdbfile)
		N.writeNative(pdbfile,CG_pos,grofile)
	if args.aa_pdb:
		#Generating sample DNA-RNA files
		#Based on the user input parameters the all_atom structure of B-form dsDNA and A-form ds-RNA (generated using NUCGEN-plus:http://nucleix.mbu.iisc.ernet.in/nucgenplus/about.html)
		# was converted into a 3 bead coarse grain form
		sample_outpdb = "sample_b-dna.pdb"
		sample_outgro = "gro_b-dna.gro"
		N.writeNative(sample_outpdb,CG_pos,sample_outgro)
		print ("EXITING")
		exit()
		pdbfile = args.aa_pdb
		if not args.all_chains:
			pdbfile = PrePDB().homoNmer2Monomer(pdbfile)
		(pdbfile,terminal_residues) = PrePDB().appendChain(pdbfile)
		N.writeNative(pdbfile,CG_pos,grofile)
	if args.grotop and args.aa_pdb :
		topfile = args.grotop
		nativefile="nuc_native_P-S-B.pdb"
		N.writeGromacsTop(topfile,pdbfile,nativefile,btparams,rad,fconst,cutoff,no_Pcharge,native_pairs,P_stretch,excl_rule)
	print ("done")
	


	#tobeCont
if __name__ == "__main__":
  main()


