#!/usr/bin/env python

"""
PyeSBM:
usage: python nucproSBM.py --options

"""
from __future__ import print_function
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO,Select
import collections
import numpy as np
import math
import mdtraj as md
from eSBM import esbm
from nucSBM import nucsbm,PrePDB
from gr import conmaps
from util import Utils

class Utility(Utils):
	def make_dir_sub_struc(self,mddir,mdfiles):
		#The function, copies and keeps the files in their respective directoriress
		import shutil
		for i in mdfiles:
			shutil.copy2(i,mddir)

class nucprosbm():
	def __init__(self):
		return
	def coordinateTransform(self,aafile,nucfile):
		#when a custom nucleic acid structure is added to the pdb,
		#the co-ordintaes of the RNA/DNA may overlap with those of protein.
		#This function moves the geometric center if the the DNA/RNA at a distance
		#equal to sum of radius of protein and DNA/RNA macromolecule
		#Radius is defined here as the distance of the farthest atom from the geometric center

		#this function is called if custom_nuc input is true
		#and returns name of the new aligned nuc pdb file
		
		#Loading aa_pdbfile and reading coordinates
		fin = open(aafile)
		coord = list()
		aalines = fin.readlines()
		for i in aalines:
			if i.startswith("ATOM"):
				coord.append(np.float_([i[30:38],i[38:46],i[46:54]]))
		fin.close()
		
		#protein coordinates
		coord = np.array(coord)
		#geometric center of aa structure
		#calculated by avereaging over all atomic co-ordibares
		aa_geocent = sum(coord)/len(coord);coord = coord-aa_geocent
		
		#maximum distance of any atom from geometric center
		#this ditacnce will define the radius of a sphere that can completely surround the protein
		aa_rad = 0
		for i in coord:
			d = sum(i**2)**0.5
			if d > aa_rad:
				aa_rad = d
		print ("Protein Geometric center = ",aa_geocent," approximate radius = ",aa_rad)

		
		#Loading nuc_pdbfile and reading coordinates
		fin = open(nucfile)
		nuclines = fin.readlines()
		fin.close()
		coord = list()
		for i in nuclines:
			if i.startswith("ATOM"):
				coord.append(np.float_([i[30:38],i[38:46],i[46:54]]))
				total_nuc_residues = int(i[22:26])
		#nucleotide coordinates
		coord = np.array(coord)
		#geometric center of aa structure
		#calculated by averaging over all coordinates
		nuc_geocent = sum(coord)/len(coord);new_coord = coord-nuc_geocent
		nuc_rad = 0
		for i in new_coord:
			d = sum(i**2)**0.5
			if d > nuc_rad:
				nuc_rad = d
		print ("DNA/RNA Geometric center = ",nuc_geocent," approximate radius = ",nuc_rad)
		
		#minimum distance to be kept between protein and DNA/RNA
		dist = aa_rad+nuc_rad
		#for diagonal distances a*(3)^0.5 = dist, hence a = dist/3**0.5
		a = dist/3**0.5
		#definfing new coordinates for nucleic acid structure
		new_nuc_geocent = aa_geocent + np.array([a,a,a])
		print ("DNA/RNA new Geometric center = ",new_nuc_geocent)
		
		#defining linear transition matrix for moving FNA/RNA
		trans_matrix = new_nuc_geocent - nuc_geocent

		#output pdb file
		outfile = "align_"+nucfile
		fout = open(outfile,"w+")
		count = 0
		#writing new co-ordinates to the file
		for i in nuclines:
			if i.startswith("ATOM"):
				details = i[:30]
				occ_and_B = i[54:]
				c = coord[count] + trans_matrix
				X = str('%8.3f' % c[0]).rjust(8);Y = str('%8.3f' % c[1]).rjust(8); Z = str('%8.3f' % c[2]).rjust(8)
				#print (i[30:54],coord[count])
				fout.write(details+X+Y+Z+occ_and_B)
				count = count + 1
			else:#TER and END
				fout.write(i)
		fout.close()
		fout = open(aafile,"w+")
		#appaending residue number in aafile
		for i in aalines:
			if not i.startswith("ATOM"):
				fout.write(i)
			else:
				resnum = str(int(i[22:26])+total_nuc_residues).rjust(len(i[22:26]))
				fout.write(i[:22]+resnum+i[26:])
		fout.close()
		return outfile
		#end of function
	def mergeNativePDB(self,aafile,nucfile):
		#This function merges the protein and the DNA/RNA native Coarde grain pdbfile
		#the function returns name of the merged native file

		#loading the two pdbfiles
		#fin1 = PrePDB().readPDB(nucfile)
		#fin2 = PrePDB().readPDB(aafile)
		file = open(nucfile); fin1 = file.readlines(); file.close()
		file = open(aafile); fin2 = file.readlines(); file.close()
		
		#output native file
		nativefile = "native_CB_P-S-B.pdb"
		fout = open(nativefile,"w+")
		#writing from file 1
		for l in fin1:
			if l.startswith("ATOM"):
				fout.write(l.strip()+"\n")
			elif l.startswith("TER"):
				#terminal
				ter = "TER".ljust(6)+l[6:11]+" ".center(5)+l[16:27]+"\n"
				fout.write(ter.strip()+"\n")
		l = fin1[len(fin1)-1]

		#atnum will be used to get the new atom numbers for the second file
		atnum = int(l[6:11].lstrip().rstrip())
		
		#writing second file
		for l in fin2:
			if l.startswith("ATOM"):
				l2 = l[0:6]+str(atnum+int(l[6:11].lstrip().rstrip())).rjust(5)+l[11:]
				fout.write(l2.strip()+"\n")
			elif l.startswith("TER"):
				ter = "TER".ljust(6)+l2[6:11]+" ".center(5)+l2[16:27]+"\n"
				fout.write(ter)
		fout.write("END")
		fout.close()
		return nativefile
		#end of function
	def nucproContacts(self,pdbfile,nativefile,cutoff,	strength):
		#calculating native cross contacts (if the DNA/RNA structure if native)
		print (">>Writing Nuc-amino_acid contacts<<")
		cutoff=float(cutoff)
		print (">>>cutoff= ",cutoff," A")
		#loading all atom pdb file
		load_atom = PDBParser().get_structure("test",pdbfile).get_atoms()
		#storing atom objects in a list
		atom = [a for a in load_atom]
		#identifying all atom contacts
		#Note: unlike protein CG models, nucleotide models should incorporate interactions between adj. residues (resi2-resi1>=1)
		nuc_res = ("A","G","T","U","C","DA","DG","DT","DC")
		aa_res  = esbm().amino_acid_dict2()
		nucpro_contacts = dict()
		natoms = len(atom)
		print ("Determining contacts")
		for i in atom:
			r1 = i.get_parent().id[1]
			rname1 = i.get_parent().get_resname().strip()
			#if residue is nucleotide
			if rname1 in nuc_res:
				for j in atom:
					d = float(j-i)#distance
					r2 = j.get_parent().id[1]
					rname2 = j.get_parent().get_resname().strip()
					#if residue is aminoacid and distance is less that cutoff
					if d <= cutoff and rname2 in aa_res:
						if ((r2,r1),(j.get_id(),i.get_id())) not in nucpro_contacts:
							nucpro_contacts[(r1,r2),(i.get_id(),j.get_id())] = [d,(rname1,rname2)]
		print ("Total nucleotide-protein interface contacs = ",len(nucpro_contacts))
		#loading nativefile (CG model file)
		CG_file = PrePDB().readPDB(nativefile)
		CG_coord = dict()
		for l in CG_file:
			key = (int(l[22:26].strip()),l[12:16].strip())
			X = float(l[30:38].rstrip().lstrip());Y = float(l[38:46].rstrip().lstrip());Z =	float(l[46:54].rstrip().lstrip())	
			A = int(l[6:11].strip())
			value = np.array([X,Y,Z])
			CG_coord[key] = [value,A]
			#print (key,value)
		#for k,v in CG_coord.items():
		#	print (k,v[0])
		#defining atoms in CG bead
		atoms_in = dict()
		atoms_in = PrePDB().atomIn()
		aromatic_contacts = dict()	#bases are stacking if stleast 6*(6-1) = 30 in the interacting ditances
		nucpro_cg = dict()
		for pair,distance in nucpro_contacts.items():
			pair_id  = pair[0]
			pair_atom = pair[1]
			pair_resname = distance[1]
			pair_distance = distance[0]
			a1 = [group for group in atoms_in if pair_atom[0] in atoms_in[group]]
			#a1 is list of 1 element Therefore
			a1 = a1[0]
			if pair_atom[1] in atoms_in["CA"]:
				a2 = "CA"
			else:
				a2 = "CB"
			if a1 == "S":
				if pair_resname[0][0] == "D":
					a1 = "D" + a1	#adding D for deoxyribose
				else:
					a1 = "R" + a1	#adding R for ribose
			elif a1 == "B":
				if pair_resname[0][0] == "D":
					a1 = a1 + pair_resname[0][1]  #Adding A/T/G/C/
				else:
					a1 = a1 + pair_resname[0][0]  #Adding A/T/G/c/u
			pair_atom = (a1,a2)
			#print (pair_atom)
			#stacking aromatic contacts
			#print (a1[0],pair_resname[1])
			if a1[0]=="B" and pair_resname[1] in ["HIS","TYR","PHE","TRP"]:
				if (pair_id,pair_atom) not in aromatic_contacts:
					aromatic_contacts[(pair_id,pair_atom)] = 0
				aromatic_contacts[(pair_id,pair_atom)] = aromatic_contacts[(pair_id,pair_atom)] + 1
			#invoking CG coordinates from nativefile
			#storing CG coordinates from nativefile and calculating distances
			coord1 = CG_coord[pair_id[0],pair_atom[0]][0]
			coord2 = CG_coord[pair_id[1],pair_atom[1]][0]
			CG_dist = (np.sum((coord2-coord1)**2))**0.5 
			#storing atom id grom the nativefile
			atom1_id = CG_coord[pair_id[0],pair_atom[0]][1]
			atom2_id = CG_coord[pair_id[1],pair_atom[1]][1]
			#print (pair_id,pair_atom,CG_dist)	
			nucpro_cg[(pair_id,pair_resname,pair_atom)] = (CG_dist,atom1_id,atom2_id)

		contacttype = 2
		#Writing aromatic contacts
		print (">>Aromatic contacts<<")
		for k,v in aromatic_contacts.items():
			print (k,v)
		#Storing a list of contacts for pair section
 		CG_contacts = list()
		for key,value in nucpro_cg.items():
			distance = value[0]; atom1_id = value[1]; atom2_id = value[2]
			a1 = key[2][0]; a2 = key[2][1]
			#atom1_id,atom2_id,contacttype,distance,epsilon, residue 1, atomname1, residue 2, atomname2
			if (key[0],key[2]) in aromatic_contacts and aromatic_contacts[(key[0],key[2])] >= 30:
				aromatic_aa = aa_res[key[1][1]]
				print ("Adding stacking contact b/w residue:",k[1][0],aromatic_aa,atom1_id,atom2_id)
				CG_contacts.append((atom1_id,atom2_id,contacttype, distance, strength[key[2][0]][aromatic_aa], key[0][0], a1, key[0][1], a2))
			CG_contacts.append((atom1_id,atom2_id,contacttype, distance, 1, key[0][0], a1, key[0][1], a2))			
		CG_contacts.sort()
		#writing 
		if len(nucpro_contacts)>0:
			fout = open("nucpro_inter_contacts.t","w+")
			for pair,distance in nucpro_contacts.items():
				pair_id  = pair[0]
				pair_atom = pair[1]
				pair_resname = distance[1]
				pair_distance = distance[0]
				fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(pair_atom[0],pair_atom[1],str(pair_id[0]),str(pair_id[1]),pair_resname[0],pair_resname[1],str(pair_distance)))
			fout.close()
		return CG_contacts
	def nucproInterfaceQ(self,pdbfile,nativefile,cutoff):
		print (">>Writing Nuc-amino_acid charged contacts<<")
		#loading all atom pdb file
		load_residues = PDBParser().get_structure("test",pdbfile).get_residues()
		#storing atom objects in a list
		resi = [a for a in load_residues]
		nuc_res = ["DA","DG","DC","DT","DU","A","G","C","T","U"]
		charged_AA = ["LYS","ARG","HIS","ASP","GLU"]
		Q_AA = [r for r in resi if r.get_resname() in charged_AA ]
		Q_nuc = [r for r in resi if r.get_resname().strip() in nuc_res]

		Q_pairs = list()
		count = 0
		for aa_residue in Q_AA:
			for aa_atom in aa_residue :
				for nuc_residue in Q_nuc: 
					for nuc_atom in nuc_residue:
						raa  = aa_atom.get_parent().id[1]
						rnuc = nuc_atom.get_parent().id[1]
						if [rnuc,raa] not in Q_pairs: 
							Q_pairs.append([rnuc,raa])

		fnat = open(nativefile)
		#native atom number fo residues in pair
		native_XP_or_CB = dict()
		for i in fnat:
			if i[12:16].strip() in ["CA","XP"]:
				native_XP_or_CB[int(i[22:26])] = int(i[6:11])
		fnat.close()

		temp = list()
		for i in Q_pairs:
			try:
				temp.append([native_XP_or_CB[i[0]],native_XP_or_CB[i[1]]])
			except:
				print ("Seems like Residue",i[0],"or",i[1],"lack CB or XP")
		Q_pairs = temp
		return Q_pairs
	def interfaceParams(self,defaults,CBradii,aalist,nuc_atoms):
		#opening user defined input file for interface
		#this is used to input user defined parametes for protein-DNA/RNA interfacce

		#loading input file
		#bash generate_NucproInterface_inputfile.sh
		#for generating the input file before simulations
		fin = open("nucpro_interface.input")
		all = [x for x in fin.readlines() if len(x.strip())>0 and not x.startswith(";")]
		fin.close()

		print ("Number of paramters: ",len(all))

		#defining default parameters
		defaults["CB_rad"] = 0 #buffer value
		defaults["LJ_atr"] = False; defaults["natCharge"] = False; defaults["Aro_cont"] = True

		# list of parameters in file
		#the list defined is in specific order and so is file
		#changing the order in file or in this list can lead to termination of code or wrong output
		params_list = [ "LJ_atr", "Aro_cont",  "A" , "G", "T" , "C" , "U" , "S" , "P" , "CA_rad" , "CB_rad" , "natCharge"] 

		#checking for input
		params = dict()
		for i in range(0,len(all)):
			try:
				#input can be True/False or a dfloat
				val = (all[i].split("]")[1].strip())
				if val.lower() == "true":
					val = True
				elif val.lower() == "false":
					val = False
				else:
					val = float(val)
				#in case of invalid input type or no input exception is raised
				params[params_list[i]] = val
				print ("Using",params_list[i],"as",val)
			except:
				#using default paramters if exception is raised
				params[params_list[i]] = defaults[params_list[i]]
				print ("Using default condition for",params_list[i])
	

		#aquiring CB radius
		if params["CB_rad"] == 0:
			#i.e if user have not been defined any CB redii, use default
			if not CBradii:
				#use from the prodefined dict
				aa_rad = esbm().amino_acid_radius_dict()
			elif CBradii:
				#read from file defined earlier in the program
				aa_rad = dict()
				fin = open("radii.dat")
				for i in fin:
					f = i.split()
					aa_rad[f[0].strip()] = float(f[1])
				fin.close
	
		#Getting CB atom types as mentioned in Amino acid top file
		all_CBrads = dict()		
		#this will include aromatic static contacts if they allready exists in nonbond_params
		if params["Aro_cont"]:
			stacking_term_exists = True
			#Aromatic contacts are defined by default if aromatic contacts are true
			#hence will exclude them
			residues_to_exclude = ["HIS","PHE","TRP","TYR"]
		else:
			stacking_term_exists = False 
			residues_to_exclude = ["NILL"] #buffer term
		CB_atoms_to_exclude = list()

		for i in aalist:
			a = i.split()
			if a[4].strip() == "CB": #atom type in atoms section of topfile
				if params["CB_rad"] == 0:
					#no user input for CB_rad, using default values
					all_CBrads[a[1].strip()] = aa_rad[a[3].strip()]
				elif params["CB_rad"] != 0:
					#using user defined values
					all_CBrads[a[1].strip()] = params["CB_rad"]
			if a[3].strip() in residues_to_exclude:
				#defining list of atom_types that to be excluded
				CB_atoms_to_exclude.append(a[1].strip())
		#generating interface nonbonded parms list
		epsilonij = 1; function = 1
		nonbond_params = list()
		all_CBrads["CA"] = params["CA_rad"] #adding CA_rad for easy computation

		#loading all non_bond paramters for interface	
		for y in all_CBrads:
			try:
				sorter = int(y[2:])	#flag to sort the nonbond_params list
			except:
				sorter = 0
			
			Cj12 = ((all_CBrads[y]*0.1)**12)*5*epsilonij
			for x in params:
				Ci12 =  ((params[x]*0.1)**12)*5*epsilonij
				if "B"+x in nuc_atoms and y not in CB_atoms_to_exclude:	#for all present N bases
									#Combination rule C12 = (C12i*C12j)^(1/2)
					A = 0; B = (Ci12*Cj12)**(0.5)
					nonbond_params.append([sorter,y,"B"+x,function,A,B])
				elif "D"+x in nuc_atoms:	#for all presnt deoxy sugars represented by DS
					#Combination rule C12 = (C12i*C12j)^(1/2)
					A = 0; B = (Ci12*Cj12)**(0.5)
					nonbond_params.append([sorter,y,"D"+x,function,A,B])
				elif "X"+x in nuc_atoms:	#for phosphate represented by XP
					#Combination rule C12 = (C12i*C12j)^(1/2)
					A = 0; B = (Ci12*Cj12)**(0.5)
					nonbond_params.append([sorter,y,"X"+x,function,A,B])
				elif "R"+x in nuc_atoms:		#for ribose sugar represented by RS
					#Combination rule C12 = (C12i*C12j)^(1/2)
					A = 0; B = (Ci12*Cj12)**(0.5)
					nonbond_params.append([sorter,y,"R"+x,function,A,B])
		nonbond_params.sort()
		print ("Inter Protein-NucleicAcid nonbonded parameters")
		for i in nonbond_params:
			print (i)
		return params,nonbond_params
		#end of function
	def mergeTopfile(self,aa_topfile,nuc_topfile,cutoff,pdbfile,nativefile,rad,CBradii,interface,custom_nuc):
		#The function reads and merges AA and nucleotide topfiles
		#This requires merging section wise and accordingly changing atom numbers
		#the function returns merged topfile

		#reading stacking epsilon
		#stacking terms are represented by LJ C10-C12 terms
		stackeps = PrePDB().stackEps()
		
		#reading nuc and aa topfile
		fin1 = open(nuc_topfile)
		fin2 = open(aa_topfile)
		top1 = fin1.read()
		top2 = fin2.read()
		fin1.close(); fin2.close()

		# all headers are enclosed in []. 
		#spliting using "[" splits all section into a list 
		# headers ] data_in_the_section 
		top1 = top1.split("[")  #split sections at the start of headers
		top2 = top2.split("[")  #split sections at the start of headers
		nuclist = list(); aalist = list()

		#splitting topfile based on headers
		# headers ] data_in_the_section, hence spliting using delim = "]"
		for i in top1:
			i=i.rstrip().lstrip()
			a = i.split("]")
			if len(a)==2:#if the section exists
				nuclist.append([a[0].strip(),a[1].strip()])
		for i in top2:
			i=i.rstrip().lstrip()
			a = i.split("]")
			if len(a)==2:#if the section exists
				aalist.append([a[0].strip(),a[1].strip()])

		#writing topdile
		topfile = aa_topfile.split("_")[1]
		fout = open(topfile,"w+")
		#header
		fout.write(top2[0])
		order = [""]
		print ("Native CG Nucleotide atoms:",len(nuclist))
		print ("Native CG Amino Acid atoms:",len(aalist))
		
		#getting nucleotide atoms
		for n in nuclist:
			header = n[0]
			if header == "atomtypes":
				nuc_atoms_section = [x for x in n[1].split("\n") if not x.startswith(";")]
		nuc_atoms = list()
		for i in nuc_atoms_section:
			nuc_atoms.append(i.split()[0].strip())

		#getting Aromactic CBs
		for a in aalist:
			header = a[0]
			if header == "atoms":
				CB_list  = [x for x in a[1].split("\n") if not x.startswith(";")]	#loading only atom section in aalist
				#clearing sections that are not needed
		#by default aromatic_contacts = True
		aromatic_contacts = True
		
		if interface:
			params,nonbond_params = self.interfaceParams(rad,CBradii,CB_list,nuc_atoms)
			aromatic_contacts = params["Aro_cont"]

		aromatic_CB_Base = list()
		if aromatic_contacts:
			#defining AA-nucleotide stacking poteintiial in form of a LJ C10-C12 term
			function = 1
			aa_res  = esbm().amino_acid_dict2()
			for i in CB_list:
				a = i.split()
				y = a[1].strip() #atomname
				if a[4] == "CB" and a[3].strip() in ["TRP","TYR","PHE","HIS"]:
					for x in nuc_atoms:
						if x.startswith("B"):
							epsilonij = float(stackeps[x][aa_res[a[3].strip()]])
							#print (rad["stack"],epsilonij)
							B = ((rad["stack"]*0.1)**12)*5*epsilonij
							A = ((rad["stack"]*0.1)**10)*6*epsilonij
							aromatic_CB_Base.append([int(y[2:]),x,y,function,A,B])

		if custom_nuc:
			#if custom nucleotide is used, no native contacts will be used
			#defining non-native stacking poteintial for the base stacking 
			function = 1
			for yi in range(0,len(nuc_atoms)):
				y = nuc_atoms[yi]
				if y.startswith("B"):
					for xi in range(yi,len(nuc_atoms)):
						x = nuc_atoms[xi]
						if x.startswith("B"):
							epsilonij = float(stackeps[x][y])
							#print (rad["stack"],epsilonij)
							B = ((rad["stack"]*0.1)**12)*5*epsilonij
							A = ((rad["stack"]*0.1)**10)*6*epsilonij
							aromatic_CB_Base.append([int(0),x,y,function,A,B])

		print ("Aromatic stacking terms")
		for i in aromatic_CB_Base:
			print (i)

		for i in range(0,len(aalist)):
			header1 = nuclist[i][0]
			header2 =  aalist[i][0]
			data1 = nuclist[i][1]
			data2 =  aalist[i][1]
			if header1==header2:
				fout.write('\n%s%s%s\n' % ("[ ",header2," ]"))
				if header2 == "defaults":
					fout.write(data1+"\n")
				elif header2 == "atomtypes":
					data = data1.split("\n")
					for d in data:
						fout.write('  %s\n' % (d.strip()))
					data = data2.split("\n")
					for d in data:
						if not d.lstrip().startswith(";"):
							fout.write('  %s\n' % (d.strip()))
					#for non-bonded params:
					header2 = "nonbond_params"
					fout.write('\n%s%s%s\n' % ("[ ",header2," ]"))
					fout.write('; i\tj\tfunc\tA(C10)\tB(C12)\n')
					#use interface if parameters are different for A-A, N-N are not same as A-N
					if interface:
						for i in nonbond_params:
							fout.write('  %s\t%s\t%d\t%e\t%e\n' % (i[1],i[2],i[3],i[4],i[5]))
					fout.write("; amino acids + nucleic acis aromatic pi-pi stacking terms\n")
					for i in aromatic_CB_Base:
						fout.write('  %s\t%s\t%d\t%e\t%e\n' % (i[1],i[2],i[3],i[4],i[5]))
				elif header2 in ["moleculetype","system","molecules"]:
					fout.write(data1+"\n")
				elif header2 == "atoms":
					data = data1.split("\n")
					for d in data:
						atnum = d.strip().split()[0]
						fout.write('  %s\n' % (d.strip()))
					atnum = int(atnum.strip()) #store the last atom number as atnum
					data = data2.split("\n")
					for d in data:
						if not d.lstrip().startswith(";"):
							d = d.strip().split()
							atnum2 = int(d[0].lstrip().rstrip())
							fout.write('  %s\t%s\t%s\t%s\t%s\t%s\t%6.2f\t%6.2f\n' % (str(atnum+atnum2),d[1],d[2],d[3],d[4],str(atnum+atnum2),float(d[6]),float(d[7])))
				#since atoms section comes before bond/angle/dihedral and pairs, the atnum value stored above can be used in other section
				elif header2 in ["bonds","pairs"]:
					data = data1.split("\n")
					if header2=="bonds" or (header2=="pairs" and not custom_nuc):
						#write nucleotide part if its atom section or 
						#its pair section but cutom_nuc is not true, It doesn't make sense 
						# to include native contacts in a custom made structure, 
						for d in data:   
							fout.write('  %s\n' % (d.strip()))
					data = data2.split("\n")
					for d in data:
						if not d.lstrip().startswith(";"):
							d = d.strip().split()
							atnum2 = [int(d[0])+atnum,int(d[1])+atnum]
							fout.write('  %s\t%s\t%s\t%s\t%s\n' % (str(atnum2[0]),str(atnum2[1]),d[2],d[3],d[4]))
							function = d[2]
					#if header is [ pairs ], determine cross contact terms
					if header2=="pairs" and interface and params["LJ_atr"]:
						#if interface input parameters are given and params["LJ_atr"] are true
						print (">>Writing Nuc-amino_acid contacts<<")
						CG_contacts = self.nucproContacts(pdbfile,nativefile,cutoff,stackeps)
						for i in CG_contacts:
							epsilonij=float(i[4])
							B=((i[3]*0.1)**12)*5*epsilonij
							A=((i[3]*0.1)**10)*6*epsilonij
							fout.write('  %d\t%d\t%s\t%e\t%e\n' % (i[0],i[1],function,A,B))
					if header2=="pairs" and interface and params["natCharge"]:
						#if interface input parameters are given and params["natCharge"] are true
						Q_contacts = self.nucproInterfaceQ(pdbfile,nativefile,cutoff)		
						for i in Q_contacts:
							B=0.0
							A=0.0
							fout.write('  %d\t%d\t%s\t%e\t%e\n' % (i[0],i[1],function,A,B))
				elif header2 == "angles":
					data = data1.split("\n")
					for d in data:   
						fout.write('  %s\n' % (d.strip()))
					data = data2.split("\n")
					for d in data:
						if not d.lstrip().startswith(";"):
							d = d.strip().split()
							atnum2 = [int(d[0])+atnum,int(d[1])+atnum,int(d[2])+atnum]
							fout.write('  %s\t%s\t%s\t%s\t%s\t%s\n' % (str(atnum2[0]),str(atnum2[1]),str(atnum2[2]),d[3],d[4],d[5]))
				elif header2 == "dihedrals":
					data = data1.split("\n")
					for d in data:   
						fout.write('  %s\n' % (d.strip()))
					data = data2.split("\n")
					for d in data:
						if not d.lstrip().startswith(";"):
							d = d.strip().split()
							atnum2 = [int(d[0])+atnum,int(d[1])+atnum,int(d[2])+atnum,int(d[3])+atnum]
							fout.write('  %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (str(atnum2[0]),str(atnum2[1]),str(atnum2[2]),str(atnum2[3]),d[4],d[5],d[6],d[7]))	
				elif header2 == "exclusions":
					if not custom_nuc:	#if not using custom nucleotide structure
						data = data1.split("\n")
						for d in data:   
							fout.write('  %s\n' % (d.strip()))
					data = data2.split("\n")
					for d in data:
						if not d.lstrip().startswith(";"):
							d = d.strip().split()
							atnum2 = [int(d[0])+atnum,int(d[1])+atnum]
							fout.write('  %s\t%s\n' % (str(atnum2[0]),str(atnum2[1])))
					if interface  and params["LJ_atr"]:
						for d in CG_contacts:
							fout.write('  %s\t%s\n' % (str(d[0]),str(d[1])))
		##########Loop end######
		fout.close()
		return topfile		
	def mergeGrofile(self,aa_grofile,nuc_grofile):
		#this function merges nuc and AA .gro files
		aa_grofile = "aa_gromacs.gro"
		nuc_grofile = "nuc_gromacs.gro"
		#reading grofile
		fin1 = open(nuc_grofile)
		fin2 = open(aa_grofile)
		file1 = fin1.read().split("\n")
		file2 = fin2.read().split("\n")
		fin1.close(); fin2.close()
		#writing combined file
		grofile = aa_grofile.split("_")[1]
		fout = open(grofile,"w+")
		fout.write('%s\n' % file1[0])
		fout.write('%s\n' % str(int(file1[1].lstrip().rstrip())+int(file2[1].lstrip().rstrip())).rjust(5))
		for i in range(2,len(file1)-1):
			fout.write('%s\n' % file1[i])
		count = 0
		for i in range(2,len(file2)-1):
			count = count + 1
			fout.write('%s%s%s\n' % (file2[i][:15],str(count+int(file1[1].lstrip().rstrip())).rjust(5),file2[i][20:]))
		fout.write('%s\n' % file2[len(file2)-1])
		fout.close()
		return grofile
	def genContactfile(self,topfile,nativefile):
		print (">>Generating SMOG-like contact file")
		ftop = open(topfile)
		file = ftop.read()
		file = file.split("[")
		ftop.read
		for i in file:
			i=i.strip()
			if i.startswith("exclusions"):
				header = i.split("]")[0].strip()
				data = i.split("]")[1].split("\n")
		fnat = open(nativefile)
		chain = dict()
		prev_chain = ""
		chain_count = 0
		for l in fnat:
			if not l.startswith("END"):
				if l[21].strip() != prev_chain:
					chain_count = chain_count + 1
				#chain[int(l[6:11])] = ord(l[21])-64
				chain[int(l[6:11])] = chain_count
				prev_chain = l[21].strip()
		fnat.close()
		fout = open("Nuc-AA_contact.txt","w+")
		for i in data:
			i = i.strip()
			if not i.startswith(";") and len(i)!=0:
				a1 = int(i.split()[0].strip())
				a2 = int(i.split()[1].strip())
				c1 = chain[a1]
				c2 = chain[a2]
				fout.write('%d %d %d %d\n' % (c1,a1,c2,a2))
		return "Nuc-AA_contact.txt"
def main():
	import argparse
	parser = argparse.ArgumentParser(description="Generate GROMACS and OPTIM potential files for enhanced SBM models.")
	parser.add_argument("--CA_rad","-CA_rad", help="Radius for C-alpha atom. Default=4.0")
	parser.add_argument("--CAcom","-CAcom",action='store_true',help="Place C-alpha at COM of backbone")
	parser.add_argument("--CB_rad","-CB_rad", help="Statistically derived values by default for attype 2.")
	parser.add_argument("--Kb","-Kb", help="Kbond")
	parser.add_argument("--Ka","-Ka", help="Kangle")
	parser.add_argument("--Kd","-Kd", help="Kdihedral")
	parser.add_argument("--cutoff","-cutoff", help="Cut-off for contact-map generation")
	parser.add_argument("--scaling","-scaling", help="Scaling for mapping to all-atom contact-map.")
	parser.add_argument("--attype", "-attype",help="Number of atom types. E.g. 1 for CA, 2 for CA and CB")
	parser.add_argument("--CBcom","-CBcom", action='store_true', default=False,help='Put CB at center of mass of side-chain (no hydrogens)')
	parser.add_argument("--CBfar", "-CBfar", action='store_true', help="Place C-beta on farthest non-hydrogen atom.")
	parser.add_argument("--dsb", "-dsb",action='store_true', help="Use desolvation barrier potential for contacts.")
	parser.add_argument("--native_ca","-native_ca", help='Native file with only C-alphas. Just grep pdb. ')
	parser.add_argument("--grotop","-grotop",help='Gromacs topology file output name.')
	parser.add_argument("--aa_pdb","-aa_pdb", help='all-atom pdbfile e.g. 1qys.pdb')
	parser.add_argument("--pdbgro","-pdbgro", help='Name for .gro file.')
	parser.add_argument("--w_native","-w_native", help='Write native files, CA-CB_P-S-B from all atom PDB file.')
	parser.add_argument("--pl_map","-pl_map", action='store_true', default=False, help='Plot contact map for two bead model')
	parser.add_argument("--skip_glycine","-skip_glycine", action='store_true', default=False, help='Skip putting Cbeta on glycine')
	parser.add_argument('--btmap',"-btmap", action='store_true', help='Use Betancourt-Thirumalai interaction matrix.')
	parser.add_argument('--mjmap',"-mjmap", action='store_true', help='Use Miyazawa-Jernighan interaction matrix.')
	parser.add_argument('--CB_radii',"-CB_radii",action='store_true', help='External contact map in format chain res chain res')

	parser.add_argument("--P_rad", help="Radius for Backbone Phosphate group bead. Default=3.7A")
	parser.add_argument("--S_rad", help="Radius for Backbone Sugar group bead. Default=3.7A")
	parser.add_argument("--Bpu_rad", help="Radius for N-Base Purine bead. Default=1.5A")
	parser.add_argument("--Bpy_rad", help="Radius for N-Base Pyrimidine bead. Default=1.5A")
	parser.add_argument("--nKb", help="Kbond for RNA/DNA")
	parser.add_argument("--nKa", help="Kangle for RNA/DNA. Default=20")
	parser.add_argument("--nKd", help="Kdihedral for Bi-Si-Si+1-Bi+1. Default=0.5")
	parser.add_argument("--P_nKd", help="Kdihedral for Backbone Pi-Pi+1-Pi+2-Pi+3. Default=0.7")
	
	parser.add_argument("--Bpu_pos", help="Put input atom of Purine [N1,C2,H2-N2,N3,C4,C5,C6,O6-N6,N7,C8,N9,COM] as position of B. Default=COM(Center_of_Mass)")
	parser.add_argument("--Bpy_pos", help="Put input atom of Pyrimidine [N1,C2,O2,N3,C4,O4-N4,C5,C6,COM] as position of B. Default=COM(Center_of_Mass)")
	parser.add_argument("--S_pos", help="Put input atom of Sugar [C1',C2',C3',C4',C5',H2'-O2',O3',O4',O5',COM] as position of S   . Default=COM(Center_of_Mass)")
	parser.add_argument("--P_pos", help="Put input atom of Phosphate [P,OP1,OP2,O5',COM] group as position of P. Default=COM(Center_of_Mass)")
	
	parser.add_argument("--all_chains", action='store_true', default=False, help='Will Not remove identical chains.')	
	parser.add_argument("--pistacklen", help="pi-pi stacking length. Default=3.6A")

	parser.add_argument("--debye",action='store_true', help="Use debye electrostatic term")
	parser.add_argument("--T", help="System temperature. Default = 100K")
	parser.add_argument("--CBcharge","-CBcharge", action='store_true', default=False, help='Put charges on CB for K,L,H,D,E')
	parser.add_argument("--no_Pcharge","-no_Pcharge", action='store_true', default=False, help='No negative charge on Phosphate bead')
	parser.add_argument("--Kr", help="Krepulsion. Default=5.7A")
	parser.add_argument("--iconc", help="Solvant ion conc. Default=0.1M")  
	parser.add_argument("--irad", help="Solvant ion rad. Default=1.4A")  
	parser.add_argument("--dielec", help="Dielectric constant of solvant. Default=70")
	
	#disabled for now
	parser.add_argument('--hpstrength',"-hpstrength",help='Strength with which hydrophobic contacts interact.')
	parser.add_argument('--ext_conmap',"-ext_conmap",help='External contact map in format chain res chain res')
	parser.add_argument("--interaction","-interaction",action='store_true', default=False, help='User defined interactions in file interaction.dat.')
	parser.add_argument("--dswap","-dswap", action='store_true', default=False, help='For domain swapping runs. Symmetrised SBM is generated.')
	parser.add_argument('--hphobic',"-hphobic",action='store_true',help='Generate hydrophobic contacts.')
	parser.add_argument('--hpdist', "-hpdist", help='Equilibrium distance for hydrophobic contacts.')

	parser.add_argument("--interface","-interface", action='store_true', default=False, help='Takes input for Nucleiotide_Protein interface from file nucpro_interface.input.')
	parser.add_argument("--custom_nuc", help='Use custom non native DNA/RNA structure1.')


	args = parser.parse_args()


	X = esbm()
	Y = conmaps()
	U = Utils()
	N = nucsbm()
	Z = nucprosbm()
	import json
	d=X.amino_acid_radius_dict()
	with open('file.txt', 'w') as file:
		file.write(json.dumps(d))
	#Set default parameters for proteins
	atomtypes = 2
	interface = False
	custom_nuc = False
	sopc=True
	btparams=False
	mjmap=False
	btmap=False
	dswap=False
	skip_glycine=False
	Kb = 200
	Ka = 40
	Kd = 1
	cutoff = 4.5
	CA_rad = 4.0
	CAcom = False
	hphobic=False
	dsb=False
	CBfar=False
	CBcom=True
	CBcharge = False
	no_Pcharge = False
	#Set default parameters for nucleotides
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
	iconc = 0.1						#M
	D =	70							#dilectric constant (dimensionlsess)
	debye = False
	T = 120
	pur_atom = ("N1","C2","H2-N2","N3","C4","C5","C6","O6-N6","N7","C8","N9","COM")
	pyr_atom = ("N1","C2","O2","N3","C4","O4-N4","C5","C6","H7-C7","COM")
	sug_atom = ("C1'","C2'","C3'","C4'","C5'","H2'-O2'","O3'","O4'","O5'","COM")
	phos_atom = ("P","OP1","OP2","O5'","COM")
	Bpu_pos = "COM"		#Center of Mass for purine
	Bpy_pos = "COM"		#Center of Mass for pyrimidine
	P_pos = "COM"			#Center of Mass for phosphate group
	S_pos = "COM"			#Center of Mass for sugar
	#setting common parameters
	#replacing default parameters with input parameters for proteins
	if args.custom_nuc:
		custom_nuc = True
		custom_nuc_file = args.custom_nuc
	if args.interface:
		interface = True
	if args.Kb:
		Kb=float(args.Kb)
	if args.Kd:
		Kd=float(args.Kd)
	if args.Ka:
		Ka=float(args.Ka)
	if args.cutoff:
		cutoff=float(args.cutoff)
	if args.interaction:
		U.file_exists('interaction.dat')
		btparams=True
	else:
		btparams=False
	if args.CB_radii:
		U.file_exists('radii.dat')
		CBradii=True
	else:
		CBradii=False
	if args.skip_glycine:
		skip_glycine=True
		sopc = False
		#assert(args.attype==2)
	if args.CBcharge:
		CBcharge = True
	if args.no_Pcharge:
		no_Pcharge = True
	if args.CAcom:
		CAcom=True
	if args.dswap:
		dswap=True
		print ('This one assumes both chains are identical.')
	if args.CA_rad:
		CA_rad=float(args.CA_rad)
		print ("Setting CA_rad to ",CA_rad)
	if args.CB_rad:
		CB_rad = float(args.CB_rad)
		aa_resi = X.amino_acid_dict2()
		fout = open("radii.dat","w+")
		print ("CB radius given via user input. Storing in radii.dat")
		for i in aa_resi:
			fout.write('%s%4.2f\n' % (i.ljust(4),CB_rad))
		fout.close()
		CBradii = True		
	if args.hphobic:
		hphobic=True
	if args.dsb:
		dsb=True
	if args.hpdist:
		hpdist=args.hpdist
	else:
		hpdist=5.5
	if args.hpstrength:
		hpstrength = args.hpstrength
	else:
		hpstrength=1
	#set global variables
	if args.CBfar:
		CBcom=False
		CBfar=True

	if args.btmap:
		import shutil
		shutil.copy2('btmap.dat','interaction.dat')
		btparams=True;btmap=True
	if args.mjmap:
		import shutil
		shutil.copy2('mjmap.dat', 'interaction.dat')
		btparams=True;mjmap=True

	#Replacing default paramteres with input paramters for nucleotide
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
	if args.debye:
		debye = True
	if args.T:
		T = float(args.T)
	sol_params = {"conc":iconc,"rad":irad,"D":D}
	#initializing global paramters

	X.globals(Ka,Kb,Kd,CA_rad,skip_glycine,sopc,dswap,btparams,CAcom,hphobic,hpstrength,hpdist,dsb,mjmap,btmap,CBfar,CBcharge)
	#########################
	#if not args.ext_conmap and args.aa_pdb:
	#	pdbfile=args.aa_pdb
	#	Y.all_atom_contacts(pdbfile,cutoff,1.2)
	print (">>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<")
	#defining grofile
	if args.pdbgro:
		grofile = str(args.pdbgro)
	else:
		grofile = "gromacs.gro"
	nuc_grofile = "nuc_"+grofile
	aa_grofile = "aa_"+grofile
	#writing native file
	tablefile = N.writeTablefile(debye,D,iconc,irad,T)

	if args.w_native:
		pdbfile=args.w_native
		#creating separate nuc and aa pdb file
		(nuc_pdbfile,aa_pdbfile) = PrePDB().sepNucPro(pdbfile)
		Y.check_hetatm(pdbfile)
		if not skip_glycine:
			Y.check_glycineh(aa_pdbfile,False)
		#esbm
		X.write_CB_to_native(aa_pdbfile,False)
		#nucsbm
		N.writeNative(nuc_pdbfile,CG_pos,grofile)
	if args.aa_pdb:
		pdbfile=args.aa_pdb
		if not args.all_chains:
			#converting homodimer to monomer
			pdbfile = PrePDB().homoNmer2Monomer(pdbfile)
		#append chains without changing the chain id
		(pdbfile,terminal_residues) = PrePDB().appendChain(pdbfile)
		#seprating pdbs
		(nuc_pdbfile,aa_pdbfile) = PrePDB().sepNucPro(pdbfile)
		if custom_nuc:
			#realign input nuc_pdbfile to avoid overlaping coo-rdinates
			nuc_pdbfile = Z.coordinateTransform(aa_pdbfile,custom_nuc_file)
		X.write_CB_to_native(aa_pdbfile,sopc)
		if args.attype:
			atomtypes = int(args.attype)
		else:
			atomtypes = 2
			print (">>>Protein Atomtypes not defined. Using default atom types = 2 [CA,CB]")		
			#checking for input top file name
		if args.grotop:
			topfile = str(args.grotop)
		else:
			topfile = 'gromacs.top'
		nuc_topfile = "nuc_"+topfile
		aa_topfile = "aa_"+topfile
		#cheking for atom types
		if atomtypes==2:
			nativefile='native_cb.pdb'
			X.write_gromacs_top(aa_topfile,int(atomtypes),aa_pdbfile,nativefile,CA_rad,sopc,btparams,Ka,Kb,Kd,cutoff,CBcom,CBradii)
		if atomtypes==1:
			U=Utils()
			if U.file_exists('native_ca.pdb'):
				nativefile='native_ca.pdb'
				X.write_gromacs_top(topfile, int(atomtypes), pdbfile, nativefile,CA_rad,sopc,btparams,Ka,Kb,Kd,cutoff,CBcom,CBradii)
			else:
				print ('Need native_ca.pdb file. Just grep "CA" pdbfile.')
		#writing grofile for protein
		X.write_gro_gro(aa_grofile,aa_pdbfile,atomtypes,sopc)	
		if args.pl_map:
			l=np.loadtxt('contacts.txt')
			xi=l[:,1]
			yi=l[:,3]
			Y.plot_map(xi,yi,'conmap','Res2','Res1')
		#writing nucleotide gro and native file
		N.writeNative(nuc_pdbfile,CG_pos,nuc_grofile)
		nativefile="nuc_native_P-S-B.pdb"
		Stats = N.writeGromacsTop(nuc_topfile,nuc_pdbfile,nativefile,btparams,rad,fconst,cutoff,no_Pcharge)
		if Stats == 1:
			nativefile = Z.mergeNativePDB("native_cb.pdb","nuc_native_P-S-B.pdb")
			rad["CA_rad"] = CA_rad
			topfile = Z.mergeTopfile(aa_topfile,nuc_topfile,cutoff,pdbfile,nativefile,rad,CBradii,interface,custom_nuc)
			grofile = Z.mergeGrofile(aa_grofile,nuc_grofile)
			contfile= Z.genContactfile(topfile,nativefile)
		else:
			topfile = aa_topfile; grofile = aa_grofile
		U.make_dir('PATH');U.make_dir('MD')
		U.make_dir('MD/AA');U.make_dir('MD/Nuc');U.make_dir('MD/Nuc_AA')
		U.make_dir('Native_PDBs')
		#U.make_dir_struc('PATH','MD')
		lU = Utility()	#local inheritance of Utils
		lU.make_dir_sub_struc('MD/AA',[aa_grofile,aa_topfile,tablefile])
		if Stats == 1:
			lU.make_dir_sub_struc('MD/Nuc',[nuc_grofile,nuc_topfile,tablefile])
			lU.make_dir_sub_struc('Native_PDBs',["native_CB_P-S-B.pdb","nuc_native_P-S-B.pdb"])
		lU.make_dir_sub_struc('Native_PDBs',["native_cb.pdb","native_ca.pdb"])
		lU.make_dir_sub_struc('MD/Nuc_AA',[grofile,topfile,tablefile])

if __name__ == '__main__':
	main()